!=====================================================================
!  SPECSIM3D  -  Spectral simulation of 3D Gaussian random fields
!                with geometric anisotropy (nested structures, nugget)
!
!  Author : Stefano Bonduà
!           DICAM - Department of Civil, Chemical, Environmental and
!           Materials Engineering, University of Bologna
!
!  Basis  : extension to 3D of SPECSIM2 (Pardo-Iguzquiza and Dowd,
!           2003, Computers & Geosciences 29, 1013-1020), which in turn
!           builds on the Fourier integral method (Pardo-Iguzquiza and
!           Chica-Olmo, 1993, Math. Geol. 25, 177-217).
!
!  Method : 1) the covariance is sampled on a periodic grid
!              (np = power of 2 per axis) using minimum-image lags;
!           2) a 3D FFT gives the spectral density S (real, because
!              C(h)=C(-h));  amplitude = sqrt(S)
!           3) uniform random phases with Hermitian symmetry are added
!           4) an inverse 3D FFT gives the real realization.
!           Var(z) = sum(S) = C(0) (up to clipped negative values).
!
!  Differences with respect to SPECSIM2 (documented on purpose):
!    * covariance anisotropy in 3D with GSLIB angles (azimuth, dip,
!      plunge) and three ranges per structure;
!    * covariance ranges are PRACTICAL ranges as in GSLIB: spherical
!      c=1-1.5r+0.5r^3 ; exponential c=exp(-3r) ; gaussian c=exp(-3r^2),
!      with r = h/ahmax  (SPECSIM2 used exp(-h/a) and exp(-h^2/a^2));
!    * self-conjugate frequencies get a random sign (variance exactly
!      equal to S) instead of cos(phase) plus the doubling of the zero
!      frequency;
!    * optional padding (ipad=1) so that the periodic embedding is exact
!      for the requested grid; ipad=0 reproduces the legacy rule
!      np = next power of 2 >= n;
!    * double precision, dynamic memory, own radix-2 FFT, no external
!      libraries; several realizations reuse the same spectrum.
!    * nz=1 gives the 2D case.
!
!  Parameter file (list-directed, one item group per line; text after the
!  values of a line is ignored, so lines can be annotated):
!    xmin ymin zmin
!    dx dy dz
!    nx ny nz
!    nugget nstruct
!    itype sill ahmax ahmin avert azimuth dip plunge   (nstruct lines)
!    iseed nreal ipad
!    output_base_name                    (no annotation on this line)
!    nout                                (number of outputs)
!    code                                (nout lines: 1 GSLIB .out, 2 binary VTK
!                                         of every realization, 3 raw float32 .bin,
!                                         4 one XML ImageData .vti with all realizations)
!    icond                               (0 unconditional, 1 conditional)
!    [conditioning block, 9 lines, only if icond = 1: see condition_field]
!  Files that end after output_base_name are read in the legacy way (outputs and
!  conditioning chosen with --bin, --vtk, --cond on the command line).
!  itype: 1 spherical, 2 exponential, 3 gaussian.
!  Angles in degrees (x = East, y = North, z = up). The rotation matrix
!  follows the GSLIB setrot formulas, re-implemented from the published
!  equations and checked numerically (test_rotation.py): azimuth is
!  clockwise from North; a POSITIVE dip tilts the major axis towards +z;
!  plunge rotates the two minor axes about the major axis (90 deg swaps
!  ahmin and avert).  Check the dip sign against the GSLIB manual before
!  comparing with GSLIB output.
!
!  Usage : specsim3d file.par [--vario variofile.par]
!  Output (set in the parameter file): code 1 : <base>.out (GSLIB grid, nreal columns,
!          x fastest); code 2 : <base>_<r>.vtk, BINARY VTK (big-endian float32, legacy
!          format, opens directly in ParaView) for every realization; code 3 :
!          <base>.bin (raw float32, realizations one after the other, x fastest);
!          code 4 : <base>.vti, ONE VTK XML ImageData file (raw binary, little-endian)
!          with all the realizations as arrays sim_1 ... sim_n (opens directly in ParaView).
!          Legacy files (no output lines): default .out + ASCII VTK of realization 1;
!          --bin : .bin + binary VTK of realization 1; --vtk : binary VTK of all.
!  Conditioning: icond = 1 in the parameter file conditions the realizations on data
!          (kriging of the residual, super block search, optional normal score
!          transform and back transform); see the manual, section on conditional simulation.
!          Legacy: --cond condfile.par.
!  Variograms: --vario variofile.par computes the experimental variograms of
!          every realization in Fortran (all node pairs, OpenMP) along grid-aligned
!          directions and writes <base>_vario.csv with the ensemble mean, sd and
!          the model variogram. variofile.par:  ndir iperreal
!                                               sx sy sz nlag   (ndir lines)
!          lag k is k*(sx,sy,sz) grid cells; iperreal=1 adds one column per
!          realization.
!=====================================================================
module sim3d_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846d0
  integer, parameter :: maxst = 10

  type :: structure
    integer  :: itype
    real(dp) :: sill, ahmax, ahmin, avert, az, dip, plg
    real(dp) :: rot(3,3)       ! GSLIB scaled rotation matrix
    real(dp) :: axes(3,3)      ! rows: principal directions * range
  end type

  integer(8) :: rng_state = 1_8

contains

  !-------------------------------------------------------------------
  subroutine setrot(s)
    type(structure), intent(inout) :: s
    real(dp) :: alpha, beta, theta, sina, cosa, sinb, cosb, sint, cost
    real(dp) :: afac1, afac2, anis1, anis2, deg
    deg = pi/180.0d0
    if (s%az >= 0.0d0 .and. s%az < 270.0d0) then
      alpha = (90.0d0 - s%az)*deg
    else
      alpha = (450.0d0 - s%az)*deg
    end if
    beta  = -s%dip*deg
    theta =  s%plg*deg
    sina = sin(alpha); cosa = cos(alpha)
    sinb = sin(beta);  cosb = cos(beta)
    sint = sin(theta); cost = cos(theta)
    anis1 = s%ahmin/s%ahmax
    anis2 = s%avert/s%ahmax
    afac1 = 1.0d0/max(anis1, 1.0d-20)
    afac2 = 1.0d0/max(anis2, 1.0d-20)
    s%rot(1,1) = cosb*cosa
    s%rot(1,2) = cosb*sina
    s%rot(1,3) = -sinb
    s%rot(2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    s%rot(2,2) = afac1*( cost*cosa + sint*sinb*sina)
    s%rot(2,3) = afac1*( sint*cosb)
    s%rot(3,1) = afac2*( sint*sina + cost*sinb*cosa)
    s%rot(3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    s%rot(3,3) = afac2*( cost*cosb)
    ! principal directions scaled by their ranges (unscaled rows * range)
    s%axes(1,:) = s%rot(1,:)*s%ahmax
    s%axes(2,:) = s%rot(2,:)/afac1*s%ahmin
    s%axes(3,:) = s%rot(3,:)/afac2*s%avert
  end subroutine setrot

  !-------------------------------------------------------------------
  function cutoff(itype) result(c)       ! distance (in ranges) beyond which
    integer, intent(in) :: itype          ! the covariance is < ~1e-3
    real(dp) :: c
    select case (itype)
    case (1); c = 1.0d0
    case (2); c = 2.3d0
    case default; c = 1.6d0
    end select
  end function cutoff

  !-------------------------------------------------------------------
  function covariance(st, nst, nugget, dx, dy, dz) result(c)
    integer, intent(in) :: nst
    type(structure), intent(in) :: st(nst)
    real(dp), intent(in) :: nugget, dx, dy, dz
    real(dp) :: c, h, hr, rx, ry, rz
    integer :: i
    c = 0.0d0
    if (dx == 0.0d0 .and. dy == 0.0d0 .and. dz == 0.0d0) c = nugget
    do i = 1, nst
      rx = st(i)%rot(1,1)*dx + st(i)%rot(1,2)*dy + st(i)%rot(1,3)*dz
      ry = st(i)%rot(2,1)*dx + st(i)%rot(2,2)*dy + st(i)%rot(2,3)*dz
      rz = st(i)%rot(3,1)*dx + st(i)%rot(3,2)*dy + st(i)%rot(3,3)*dz
      h  = sqrt(rx*rx + ry*ry + rz*rz)
      hr = h/st(i)%ahmax
      select case (st(i)%itype)
      case (1)
        if (hr < 1.0d0) then
          c = c + st(i)%sill*(1.0d0 - 1.5d0*hr + 0.5d0*hr**3)
        end if
      case (2)
        c = c + st(i)%sill*exp(-3.0d0*hr)
      case default
        c = c + st(i)%sill*exp(-3.0d0*hr*hr)
      end select
    end do
  end function covariance

  !-------------------------------------------------------------------
  function nextpow2(n) result(m)
    integer, intent(in) :: n
    integer :: m
    m = 1
    do while (m < n)
      m = m*2
    end do
  end function nextpow2

  !-------------------------------------------------------------------
  function urand() result(u)    ! Lehmer / Park-Miller 16807 (as SPECSIM2)
    real(dp) :: u
    rng_state = mod(16807_8*rng_state, 2147483647_8)
    u = real(rng_state, dp)/2147483647.0d0
  end function urand

  !-------------------------------------------------------------------
  subroutine fft1d(x, n, tw)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: x(n)
    complex(dp), intent(in) :: tw(max(n/2,1))
    integer :: i, j, k, m, len, half, step
    complex(dp) :: t, u
    if (n < 2) return
    j = 1
    do i = 1, n-1
      if (i < j) then
        t = x(i); x(i) = x(j); x(j) = t
      end if
      m = n/2
      do while (m >= 1 .and. j > m)
        j = j - m
        m = m/2
      end do
      j = j + m
    end do
    len = 2
    do while (len <= n)
      half = len/2
      step = n/len
      do i = 1, n, len
        do k = 0, half-1
          u = x(i+k)
          t = tw(k*step+1)*x(i+k+half)
          x(i+k)      = u + t
          x(i+k+half) = u - t
        end do
      end do
      len = len*2
    end do
  end subroutine fft1d

  !-------------------------------------------------------------------
  function twiddles(n, isign) result(tw)
    integer, intent(in) :: n, isign
    complex(dp) :: tw(max(n/2,1))
    integer :: k
    do k = 0, max(n/2,1)-1
      tw(k+1) = cmplx(cos(2.0d0*pi*k/n), isign*sin(2.0d0*pi*k/n), dp)
    end do
  end function twiddles

  !-------------------------------------------------------------------
  subroutine fft3d(a, nx, ny, nz, isign)      ! no scaling in either direction
    integer, intent(in) :: nx, ny, nz, isign
    complex(dp), intent(inout) :: a(nx,ny,nz)
    complex(dp), allocatable :: buf(:), tw(:)
    integer :: i, j, k
    if (nx > 1) then
      allocate(tw(max(nx/2,1))); tw = twiddles(nx, isign)
!$omp parallel private(buf, j, k)
      allocate(buf(nx))
!$omp do collapse(2) schedule(static)
      do k = 1, nz
        do j = 1, ny
          buf = a(:,j,k); call fft1d(buf, nx, tw); a(:,j,k) = buf
        end do
      end do
!$omp end do
      deallocate(buf)
!$omp end parallel
      deallocate(tw)
    end if
    if (ny > 1) then
      allocate(tw(max(ny/2,1))); tw = twiddles(ny, isign)
!$omp parallel private(buf, i, k)
      allocate(buf(ny))
!$omp do collapse(2) schedule(static)
      do k = 1, nz
        do i = 1, nx
          buf = a(i,:,k); call fft1d(buf, ny, tw); a(i,:,k) = buf
        end do
      end do
!$omp end do
      deallocate(buf)
!$omp end parallel
      deallocate(tw)
    end if
    if (nz > 1) then
      allocate(tw(max(nz/2,1))); tw = twiddles(nz, isign)
!$omp parallel private(buf, i, j)
      allocate(buf(nz))
!$omp do collapse(2) schedule(static)
      do j = 1, ny
        do i = 1, nx
          buf = a(i,j,:); call fft1d(buf, nz, tw); a(i,j,:) = buf
        end do
      end do
!$omp end do
      deallocate(buf)
!$omp end parallel
      deallocate(tw)
    end if
  end subroutine fft3d

  !-------------------------------------------------------------------
  subroutine write_vti(fname, res, nx, ny, nz, nreal, x0, y0, z0, dx, dy, dz)
    ! VTK XML ImageData (.vti), ONE file with all realizations as point-data arrays sim_1 ... sim_n.
    ! Raw appended binary, little-endian float32, UInt64 block headers (no compression); ParaView opens it
    ! directly with the geometry (dimensions, origin, spacing) stored in the XML header.
    character(len=*), intent(in) :: fname
    integer, intent(in) :: nx, ny, nz, nreal
    real(dp), intent(in) :: res(nx,ny,nz,nreal)
    real(dp), intent(in) :: x0, y0, z0, dx, dy, dz
    character(len=1), parameter :: nl = achar(10)
    character(len=512) :: line
    character(len=40) :: ext
    real(4), allocatable :: sp(:,:,:)
    integer(8) :: nbytes, offset
    integer :: r, k, j, i
    integer(1) :: probe(4)
    ! the file is declared little-endian: stop if the machine is not
    probe = transfer(1_4, probe)
    if (probe(1) /= 1) stop 'write_vti: big-endian machines are not supported'
    nbytes = 4_8*int(nx,8)*int(ny,8)*int(nz,8)
    write(ext,'(a,i0,a,i0,a,i0)') '0 ', nx-1, ' 0 ', ny-1, ' 0 ', nz-1
    open(unit=13, file=trim(fname), status='replace', access='stream', form='unformatted')
    write(13) '<?xml version="1.0"?>'//nl
    write(13) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//nl
    write(line,'(a,a,a,3(1x,es16.8),a,3(1x,es16.8),a)') '  <ImageData WholeExtent="', trim(ext), '" Origin="', x0, y0, z0, &
          '" Spacing="', dx, dy, dz, '">'
    write(13) trim(line)//nl
    write(line,'(a,a,a)') '    <Piece Extent="', trim(ext), '">'
    write(13) trim(line)//nl
    write(13) '      <PointData Scalars="sim_1">'//nl
    do r = 1, nreal
      offset = (nbytes + 8_8)*int(r-1, 8)
      write(line,'(a,i0,a,i0,a)') '        <DataArray type="Float32" Name="sim_', r, '" format="appended" offset="', offset, '"/>'
      write(13) trim(line)//nl
    end do
    write(13) '      </PointData>'//nl
    write(13) '      <CellData/>'//nl
    write(13) '    </Piece>'//nl
    write(13) '  </ImageData>'//nl
    write(13) '  <AppendedData encoding="raw">'//nl
    write(13) '   _'
    allocate(sp(nx,ny,nz))
    do r = 1, nreal
      do k = 1, nz; do j = 1, ny; do i = 1, nx
        sp(i,j,k) = real(res(i,j,k,r), 4)
      end do; end do; end do
      write(13) nbytes
      write(13) sp
    end do
    deallocate(sp)
    write(13) nl//'  </AppendedData>'//nl
    write(13) '</VTKFile>'//nl
    close(13)
  end subroutine write_vti

  !-------------------------------------------------------------------
  subroutine write_vtk_bin(fname, a, nx, ny, nz, x0, y0, z0, dx, dy, dz, title)
    ! Legacy VTK, BINARY (big-endian float32), STRUCTURED_POINTS. x runs fastest.
    ! Uses the (non-standard but widespread) CONVERT= specifier of OPEN,
    ! supported by gfortran and Intel Fortran.
    character(len=*), intent(in) :: fname, title
    integer, intent(in) :: nx, ny, nz
    real(4), intent(in) :: a(nx,ny,nz)
    real(dp), intent(in) :: x0, y0, z0, dx, dy, dz
    character(len=1), parameter :: nl = achar(10)
    character(len=256) :: line
    open(unit=12, file=trim(fname), status='replace', access='stream', &
         form='unformatted', convert='big_endian')
    write(12) '# vtk DataFile Version 3.0'//nl
    write(12) trim(title)//nl
    write(12) 'BINARY'//nl
    write(12) 'DATASET STRUCTURED_POINTS'//nl
    write(line,'(a,3(1x,i0))') 'DIMENSIONS', nx, ny, nz
    write(12) trim(line)//nl
    write(line,'(a,3(1x,es16.8))') 'ORIGIN', x0, y0, z0
    write(12) trim(line)//nl
    write(line,'(a,3(1x,es16.8))') 'SPACING', dx, dy, dz
    write(12) trim(line)//nl
    write(line,'(a,i0)') 'POINT_DATA ', int(nx,8)*int(ny,8)*int(nz,8)
    write(12) trim(line)//nl
    write(12) 'SCALARS sim float 1'//nl
    write(12) 'LOOKUP_TABLE default'//nl
    write(12) a
    write(12) nl
    close(12)
  end subroutine write_vtk_bin

  !-------------------------------------------------------------------
  function npair_count(nx, ny, nz, lx, ly, lz) result(np)
    integer, intent(in) :: nx, ny, nz, lx, ly, lz
    integer(8) :: np
    if (abs(lx) >= nx .or. abs(ly) >= ny .or. abs(lz) >= nz) then
      np = 0
    else
      np = int(nx-abs(lx), 8)*int(ny-abs(ly), 8)*int(nz-abs(lz), 8)
    end if
  end function npair_count

  !-------------------------------------------------------------------
  subroutine vario_lag(f, nx, ny, nz, lx, ly, lz, g)
    ! semivariogram of the regular-grid field f for the integer lag vector
    ! (lx,ly,lz), using ALL node pairs: g = 0.5*mean((f(x+l)-f(x))**2)
    integer, intent(in) :: nx, ny, nz, lx, ly, lz
    real(dp), intent(in) :: f(nx,ny,nz)
    real(dp), intent(out) :: g
    integer :: i, j, k, i0, i1, j0, j1, k0, k1
    real(dp) :: acc, d
    integer(8) :: np
    np = npair_count(nx, ny, nz, lx, ly, lz)
    if (np == 0) then
      g = 0.0d0; return
    end if
    i0 = max(1, 1-lx); i1 = min(nx, nx-lx)
    j0 = max(1, 1-ly); j1 = min(ny, ny-ly)
    k0 = max(1, 1-lz); k1 = min(nz, nz-lz)
    acc = 0.0d0
    do k = k0, k1
      do j = j0, j1
        do i = i0, i1
          d = f(i+lx,j+ly,k+lz) - f(i,j,k)
          acc = acc + d*d
        end do
      end do
    end do
    g = 0.5d0*acc/real(np, dp)
  end subroutine vario_lag

  !-------------------------------------------------------------------
  subroutine vario_report(varfile, outbase, res, nx, ny, nz, nreal, dx, dy, dz, st, nst, nugget)
    ! Experimental variograms of every realization along the grid-aligned
    ! directions listed in varfile, ensemble mean/sd, and the model variogram.
    ! varfile:  ndir  iperreal
    !           sx sy sz nlag          (ndir lines; lag k = k*(sx,sy,sz) cells)
    character(len=*), intent(in) :: varfile, outbase
    integer, intent(in) :: nx, ny, nz, nreal, nst
    real(dp), intent(in) :: res(nx,ny,nz,nreal), dx, dy, dz, nugget
    type(structure), intent(in) :: st(nst)
    integer :: ndir, iper, id, il, r, nlagmax, t, ntask, ios
    integer, allocatable :: sx(:), sy(:), sz(:), nlag(:), tid(:), tlag(:)
    real(dp), allocatable :: gam(:,:,:)
    real(dp) :: c0, gm, gsd, gth, h, rms, mx, e, hx, hy, hz
    integer :: nfit
    integer(8) :: np
    character(len=512) :: fname
    integer :: i

    open(unit=20, file=trim(varfile), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'cannot open variogram file: ', trim(varfile); return
    end if
    read(20,*) ndir, iper
    allocate(sx(ndir), sy(ndir), sz(ndir), nlag(ndir))
    do id = 1, ndir
      read(20,*) sx(id), sy(id), sz(id), nlag(id)
    end do
    close(20)
    nlagmax = maxval(nlag)
    ntask = sum(nlag)
    allocate(tid(ntask), tlag(ntask), gam(nlagmax, ndir, nreal))
    gam = 0.0d0
    t = 0
    do id = 1, ndir
      do il = 1, nlag(id)
        t = t + 1; tid(t) = id; tlag(t) = il
      end do
    end do

    do r = 1, nreal
!$omp parallel do schedule(dynamic) private(t, id, il)
      do t = 1, ntask
        id = tid(t); il = tlag(t)
        call vario_lag(res(:,:,:,r), nx, ny, nz, il*sx(id), il*sy(id), il*sz(id), gam(il,id,r))
      end do
!$omp end parallel do
    end do

    c0 = nugget
    do i = 1, nst
      c0 = c0 + st(i)%sill
    end do

    fname = trim(outbase)//'_vario.csv'
    open(unit=21, file=trim(fname), status='replace', action='write')
    write(21,'(a)', advance='no') 'dir,sx,sy,sz,lag,h,npairs,gamma_model,gamma_mean,gamma_sd'
    if (iper == 1) then
      do r = 1, nreal
        write(21,'(a,i0)', advance='no') ',gamma_r', r
      end do
    end if
    write(21,'(a)') ''
    write(*,'(a)') ' variograms (Fortran): direction  lags  RMS(mean-model)  max|err|   [rising part, gamma_model < 0.95 sill]'
    do id = 1, ndir
      rms = 0.0d0; mx = 0.0d0; nfit = 0
      do il = 1, nlag(id)
        np = npair_count(nx, ny, nz, il*sx(id), il*sy(id), il*sz(id))
        if (np == 0) cycle
        hx = il*sx(id)*dx; hy = il*sy(id)*dy; hz = il*sz(id)*dz
        h = sqrt(hx*hx + hy*hy + hz*hz)
        gth = c0 - covariance(st, nst, nugget, hx, hy, hz)
        gm = sum(gam(il,id,:))/real(nreal, dp)
        if (nreal > 1) then
          gsd = sqrt(sum((gam(il,id,:)-gm)**2)/real(nreal-1, dp))
        else
          gsd = 0.0d0
        end if
        write(21,'(i0,3(",",i0),",",i0,",",es14.6,",",i0,3(",",es14.6))', advance='no') &
              id, sx(id), sy(id), sz(id), il, h, np, gth, gm, gsd
        if (iper == 1) then
          do r = 1, nreal
            write(21,'(",",es14.6)', advance='no') gam(il,id,r)
          end do
        end if
        write(21,'(a)') ''
        if (gth < 0.95d0*c0) then
          e = gm - gth
          rms = rms + e*e; mx = max(mx, abs(e)); nfit = nfit + 1
        end if
      end do
      if (nfit > 0) rms = sqrt(rms/real(nfit, dp))
      write(*,'(a,i3,3i4,i6,2f12.5)') '   dir', id, sx(id), sy(id), sz(id), nlag(id), rms, mx
    end do
    close(21)
    write(*,'(a,a)') ' variogram table: ', trim(fname)
  end subroutine vario_report

end module sim3d_mod

!=====================================================================
!=====================================================================
!  Conditional simulation module (added in SPECSIM3D v1.1)
!
!  Method: kriging of the residual (conditioning by kriging).
!    z_c(x) = z*_SK(x; data) + [ z_s(x) - z*_SK(x; z_s at data nodes) ]
!  with the same simple-kriging weights for the data and for the
!  unconditional spectral realization z_s.  Data are snapped to the
!  nearest grid node (as in GSLIB sgsim); data that fall in the same
!  node are averaged.  Neighbours are found with a GSLIB-style super
!  block search; the kriging systems are solved by Cholesky
!  factorization in an OpenMP loop over the grid nodes.  Optionally the
!  data are first normal-score transformed (GSLIB nscore) and the
!  conditioned Gaussian field is back-transformed (GSLIB backtr).
!=====================================================================
module cond_mod
  use sim3d_mod
  implicit none
  real(dp), parameter :: sqrt2 = 1.4142135623730951d0
  real(dp), parameter :: sqrt2pi = 2.5066282746310002d0

contains

  !-------------------------------------------------------------------
  function gcum(x) result(p)          ! standard normal cdf
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5d0*erfc(-x/sqrt2)
  end function gcum

  !-------------------------------------------------------------------
  function gauinv(p0) result(x)       ! inverse normal cdf (Acklam + one Halley step)
    real(dp), intent(in) :: p0
    real(dp) :: x, p, q, r, e, u
    real(dp), parameter :: a(6) = [-3.969683028665376d+01, 2.209460984245205d+02, &
         -2.759285104469687d+02, 1.383577518672690d+02, -3.066479806614716d+01, 2.506628277459239d+00]
    real(dp), parameter :: b(5) = [-5.447609879822406d+01, 1.615858368580409d+02, &
         -1.556989798598866d+02, 6.680131188771972d+01, -1.328068155288572d+01]
    real(dp), parameter :: c(6) = [-7.784894002430293d-03, -3.223964580411365d-01, &
         -2.400758277161838d+00, -2.549732539343734d+00, 4.374664141464968d+00, 2.938163982698783d+00]
    real(dp), parameter :: d(4) = [7.784695709041462d-03, 3.224671290700398d-01, &
         2.445134137142996d+00, 3.754408661907416d+00]
    p = min(max(p0, 1.0d-15), 1.0d0 - 1.0d-15)
    if (p < 0.02425d0) then
      q = sqrt(-2.0d0*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0d0)
    else if (p <= 1.0d0 - 0.02425d0) then
      q = p - 0.5d0
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0d0)
    else
      q = sqrt(-2.0d0*log(1.0d0-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0d0)
    end if
    e = 0.5d0*erfc(-x/sqrt2) - p
    u = e*sqrt2pi*exp(0.5d0*x*x)
    x = x - u/(1.0d0 + 0.5d0*x*u)
  end function gauinv

  !-------------------------------------------------------------------
  subroutine heapsort_idx(key, idx, n)       ! sorts idx so that key(idx) is ascending
    integer, intent(in) :: n
    real(dp), intent(in) :: key(n)
    integer, intent(inout) :: idx(n)
    integer :: i, tmp
    do i = n/2, 1, -1
      call sift(i, n)
    end do
    do i = n, 2, -1
      tmp = idx(1); idx(1) = idx(i); idx(i) = tmp
      call sift(1, i-1)
    end do
  contains
    subroutine sift(root, last)
      integer, intent(in) :: root, last
      integer :: parent, child, v
      parent = root
      v = idx(parent)
      do
        child = 2*parent
        if (child > last) exit
        if (child < last) then
          if (key(idx(child+1)) > key(idx(child))) child = child + 1
        end if
        if (key(idx(child)) > key(v)) then
          idx(parent) = idx(child); parent = child
        else
          exit
        end if
      end do
      idx(parent) = v
    end subroutine sift
  end subroutine heapsort_idx

  !-------------------------------------------------------------------
  function powint(xlow, xhigh, ylow, yhigh, xval, pw) result(y)
    real(dp), intent(in) :: xlow, xhigh, ylow, yhigh, xval, pw
    real(dp) :: y
    if (xhigh - xlow < 1.0d-20) then
      y = 0.5d0*(ylow + yhigh)
    else
      y = ylow + (yhigh - ylow)*((xval - xlow)/(xhigh - xlow))**pw
    end if
  end function powint

  !-------------------------------------------------------------------
  subroutine nscore(nd, v, w, g, nt, vr, gt)    ! GSLIB nscore (weighted)
    integer, intent(in) :: nd
    real(dp), intent(in) :: v(nd), w(nd)
    real(dp), intent(out) :: g(nd), vr(nd), gt(nd)
    integer, intent(out) :: nt
    integer :: i, idx(nd)
    real(dp) :: sw, cp, oldcp
    idx = [(i, i = 1, nd)]
    call heapsort_idx(v, idx, nd)
    sw = sum(w)
    cp = 0.0d0; oldcp = 0.0d0
    do i = 1, nd
      cp = cp + w(idx(i))/sw
      gt(i) = gauinv(0.5d0*(cp + oldcp))
      vr(i) = v(idx(i))
      g(idx(i)) = gt(i)
      oldcp = cp
    end do
    nt = nd
  end subroutine nscore

  !-------------------------------------------------------------------
  function backtr(y, nt, vr, gt, zmin, zmax, lt, ltp, ut, utp) result(z)   ! GSLIB backtr
    integer, intent(in) :: nt, lt, ut
    real(dp), intent(in) :: y, vr(nt), gt(nt), zmin, zmax, ltp, utp
    real(dp) :: z, cdflo, cdfhi, cdfbt, lam
    integer :: lo, hi, mid
    if (y <= gt(1)) then
      cdflo = gcum(gt(1)); cdfbt = gcum(y)
      if (lt == 2) then
        z = powint(0.0d0, cdflo, zmin, vr(1), cdfbt, ltp)
      else
        z = powint(0.0d0, cdflo, zmin, vr(1), cdfbt, 1.0d0)
      end if
      if (nt == 1) z = vr(1)
    else if (y >= gt(nt)) then
      cdfhi = gcum(gt(nt)); cdfbt = gcum(y)
      select case (ut)
      case (2)
        z = powint(cdfhi, 1.0d0, vr(nt), zmax, cdfbt, utp)
      case (4)
        lam = (vr(nt)**utp)*(1.0d0 - cdfhi)
        z = (lam/max(1.0d0 - cdfbt, 1.0d-300))**(1.0d0/utp)
      case default
        z = powint(cdfhi, 1.0d0, vr(nt), zmax, cdfbt, 1.0d0)
      end select
    else
      lo = 1; hi = nt
      do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (gt(mid) <= y) then
          lo = mid
        else
          hi = mid
        end if
      end do
      z = powint(gt(lo), gt(hi), vr(lo), vr(hi), y, 1.0d0)
    end if
  end function backtr

  !-------------------------------------------------------------------
  subroutine condition_field(res, nx, ny, nz, nreal, xmin, ymin, zmin, dx, dy, dz, &
                             st, nst, nugget, cfile, outbase, cunit)
    integer, intent(in) :: nx, ny, nz, nreal, nst
    real(dp), intent(inout) :: res(nx,ny,nz,nreal)
    real(dp), intent(in) :: xmin, ymin, zmin, dx, dy, dz, nugget
    type(structure), intent(in) :: st(nst)
    character(len=*), intent(in) :: cfile, outbase
    integer, intent(in) :: cunit        ! > 0: block of the main parameter file, already open and positioned

    ! ---- parameters from the conditioning parameter file -------------
    character(len=512) :: line, datafile
    integer :: ixc, iyc, izc, ivc, iwc, itrans, ltail, utail, ndmax
    real(dp) :: tmin, tmax, skmean, zlo, zhi, ltpar, utpar
    real(dp) :: rmax, rmed, rvert, sa1, sa2, sa3
    type(structure) :: srch
    ! ---- data -------------------------------------------------------
    integer :: nvar, nraw, ios, i, j, k, r, n, nd, nt, nout, nmerged
    real(dp), allocatable :: row(:), raw(:,:), key(:), tv(:), tw(:), tg(:), vrt(:), gtt(:)
    integer, allocatable :: ord(:), ix(:), iy(:), iz(:)
    real(dp), allocatable :: dxs(:), dys(:), dzs(:), dval(:), dwt(:), dg(:), simd(:,:)
    ! ---- super blocks -----------------------------------------------
    integer :: nsbx, nsby, nsbz, nsb, mbx, mby, mbz
    real(dp) :: bsx, bsy, bsz, ex, ey, ez
    integer, allocatable :: sbstart(:), sbidx(:), sbcount(:), sbof(:)
    integer :: t0, t1, t2, tr
    integer(8) :: nsum
    integer :: nfail, nnone, nexact, iu
    real(dp) :: c0

    call system_clock(t0, tr)
    ! ---- read conditioning parameter file ---------------------------
    if (cunit > 0) then
      iu = cunit                          ! conditioning block inside the main parameter file
    else
      iu = 15                             ! separate conditioning file (legacy, --cond)
      open(iu, file=trim(cfile), status='old', action='read', iostat=ios)
      if (ios /= 0) then
        write(*,*) 'Conditioning parameter file not found: ', trim(cfile); stop 1
      end if
    end if
    read(iu,'(a)') line
    line = adjustl(line); datafile = line(1:index(line//' ', ' ')-1)
    read(iu,*) ixc, iyc, izc, ivc, iwc
    read(iu,*) tmin, tmax
    read(iu,*) itrans, skmean
    read(iu,*) zlo, zhi
    read(iu,*) ltail, ltpar
    read(iu,*) utail, utpar
    read(iu,*) rmax, rmed, rvert
    read(iu,*) sa1, sa2, sa3
    read(iu,*) ndmax
    if (cunit <= 0) close(iu)
    if (ndmax < 1 .or. ndmax > 100) stop 'cond: ndmax must be in 1..100'
    if (rmax <= 0 .or. rmed <= 0 .or. rvert <= 0) stop 'cond: search radii must be positive'
    srch%itype = 1; srch%sill = 1.0d0
    srch%ahmax = rmax; srch%ahmin = rmed; srch%avert = rvert
    srch%az = sa1; srch%dip = sa2; srch%plg = sa3
    call setrot(srch)
    ex = sqrt(sum(srch%axes(:,1)**2)); ey = sqrt(sum(srch%axes(:,2)**2)); ez = sqrt(sum(srch%axes(:,3)**2))

    ! ---- read the data (GSLIB format) -------------------------------
    open(16, file=trim(datafile), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Data file not found: ', trim(datafile); stop 1
    end if
    read(16,'(a)') line
    read(16,*) nvar
    do i = 1, nvar
      read(16,'(a)') line
    end do
    nraw = 0
    do
      read(16,'(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) nraw = nraw + 1
    end do
    rewind(16)
    read(16,'(a)') line; read(16,*) nvar
    do i = 1, nvar
      read(16,'(a)') line
    end do
    allocate(row(nvar), raw(5, nraw), key(nraw))
    n = 0; nout = 0
    do i = 1, nraw
      read(16,'(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read(line,*) row
      if (row(ivc) < tmin .or. row(ivc) >= tmax) cycle
      j = nint((row(ixc) - xmin)/dx) + 1
      k = nint((row(iyc) - ymin)/dy) + 1
      r = 1
      if (izc > 0 .and. nz > 1) r = nint((row(izc) - zmin)/dz) + 1
      if (j < 1 .or. j > nx .or. k < 1 .or. k > ny .or. r < 1 .or. r > nz) then
        nout = nout + 1; cycle
      end if
      n = n + 1
      raw(1,n) = real(j, dp); raw(2,n) = real(k, dp); raw(3,n) = real(r, dp)
      raw(4,n) = row(ivc)
      raw(5,n) = 1.0d0
      if (iwc > 0) raw(5,n) = row(iwc)
      key(n) = real(j, dp) + real(nx, dp)*(real(k-1, dp) + real(ny, dp)*real(r-1, dp))
    end do
    close(16)
    if (n == 0) stop 'cond: no conditioning data inside the grid'
    ! ---- merge data falling in the same node ---------------------------
    allocate(ord(n))
    ord = [(i, i = 1, n)]
    call heapsort_idx(key, ord, n)
    allocate(ix(n), iy(n), iz(n), dval(n), dwt(n))
    nd = 0
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (key(ord(j+1)) /= key(ord(i))) exit
        j = j + 1
      end do
      nd = nd + 1
      ix(nd) = nint(raw(1,ord(i))); iy(nd) = nint(raw(2,ord(i))); iz(nd) = nint(raw(3,ord(i)))
      dwt(nd) = sum(raw(5, ord(i:j)))
      dval(nd) = sum(raw(4, ord(i:j))*raw(5, ord(i:j)))/dwt(nd)
      i = j + 1
    end do
    nmerged = n - nd
    allocate(dxs(nd), dys(nd), dzs(nd), dg(nd))
    dxs = xmin + (ix(1:nd)-1)*dx; dys = ymin + (iy(1:nd)-1)*dy; dzs = zmin + (iz(1:nd)-1)*dz
    write(*,'(a,i8,a,i8,a,i8,a)') ' conditioning data: ', nd, ' nodes used (', nmerged, &
          ' merged, ', nout, ' outside the grid)'

    ! ---- normal score transform ---------------------------------------
    allocate(vrt(nd), gtt(nd))
    nt = nd
    if (itrans == 1) then
      call nscore(nd, dval(1:nd), dwt(1:nd), dg, nt, vrt, gtt)
      open(17, file=trim(outbase)//'_nscore.trn', status='replace', action='write')
      write(17,'(a)') 'SPECSIM3D normal score transform table (value, normal score)'
      write(17,'(i0)') 2
      write(17,'(a)') 'value'; write(17,'(a)') 'nscore'
      do i = 1, nt
        write(17,'(2es17.8)') vrt(i), gtt(i)
      end do
      close(17)
      write(*,'(a,f9.4,a,f9.4)') ' normal score transform: data mean=', sum(dg)/nd, '  variance=', &
            sum((dg-sum(dg)/nd)**2)/nd
    else
      dg = dval(1:nd) - skmean
    end if
    open(18, file=trim(outbase)//'_conddata.dat', status='replace', action='write')
    write(18,'(a)') 'SPECSIM3D conditioning data used (grid nodes)'
    write(18,'(i0)') 5
    write(18,'(a)') 'x'; write(18,'(a)') 'y'; write(18,'(a)') 'z'; write(18,'(a)') 'value'; write(18,'(a)') 'gauss'
    do i = 1, nd
      write(18,'(5es17.8)') dxs(i), dys(i), dzs(i), dval(i), dg(i)
    end do
    close(18)

    ! ---- simulated values at the data nodes ---------------------------
    allocate(simd(nd, nreal))
    do r = 1, nreal
      do i = 1, nd
        simd(i, r) = res(ix(i), iy(i), iz(i), r)
      end do
    end do

    ! ---- super block search structure ---------------------------------
    nsbx = max(1, min(100, ceiling(nx*dx/(0.5d0*ex))))
    nsby = max(1, min(100, ceiling(ny*dy/(0.5d0*ey))))
    nsbz = max(1, min(100, ceiling(nz*dz/(0.5d0*ez))))
    if (nz == 1) nsbz = 1
    do while (real(nsbx,dp)*nsby*nsbz > 2.0d6)
      nsbx = max(1, nsbx*3/4); nsby = max(1, nsby*3/4); nsbz = max(1, nsbz*3/4)
    end do
    nsb = nsbx*nsby*nsbz
    bsx = nx*dx/nsbx; bsy = ny*dy/nsby; bsz = nz*dz/nsbz
    mbx = ceiling(ex/bsx); mby = ceiling(ey/bsy); mbz = ceiling(ez/bsz)
    allocate(sbstart(nsb+1), sbcount(nsb), sbidx(nd), sbof(nd))
    sbcount = 0
    do i = 1, nd
      sbof(i) = blockof(dxs(i), dys(i), dzs(i))
      sbcount(sbof(i)) = sbcount(sbof(i)) + 1
    end do
    sbstart(1) = 1
    do i = 1, nsb
      sbstart(i+1) = sbstart(i) + sbcount(i)
    end do
    sbcount = 0
    do i = 1, nd
      j = sbof(i)
      sbidx(sbstart(j) + sbcount(j)) = i
      sbcount(j) = sbcount(j) + 1
    end do
    write(*,'(a,3i5,a,i0,a)') ' super blocks:', nsbx, nsby, nsbz, '  (search ', ndmax, ' nearest data)'
    c0 = covariance(st, nst, nugget, 0.0d0, 0.0d0, 0.0d0)

    ! ---- kriging of the residual at every node (OpenMP) ---------------
    nsum = 0; nfail = 0; nnone = 0; nexact = 0
!$omp parallel do collapse(2) schedule(dynamic, 2) default(shared) private(i, j, k) &
!$omp&   reduction(+:nsum, nfail, nnone, nexact)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call krige_node(i, j, k, nsum, nfail, nnone, nexact)
        end do
      end do
    end do
!$omp end parallel do
    call system_clock(t1)

    ! ---- back transform / restore the mean ------------------------------
    if (itrans == 1) then
!$omp parallel do collapse(2) private(i, j, k, r) schedule(static)
      do r = 1, nreal
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              res(i,j,k,r) = backtr(res(i,j,k,r), nt, vrt, gtt, zlo, zhi, ltail, ltpar, utail, utpar)
            end do
          end do
        end do
      end do
!$omp end parallel do
    else if (skmean /= 0.0d0) then
      res = res + skmean
    end if
    call system_clock(t2)
    write(*,'(a,f6.1,a,i0,a,i0,a)') ' mean neighbours per node: ', real(nsum, dp)/(real(nx,dp)*ny*nz), &
          '  (nodes without data: ', nnone, ', data nodes: ', nexact, ')'
    if (nfail > 0) write(*,'(a,i0)') ' WARNING: kriging systems regularized: ', nfail
    write(*,'(a,f8.2,a)') ' time conditioning: ', real(t1-t0)/real(tr), ' s (kriging) + ' // &
          trim(adjustl(fmt1(real(t2-t1)/real(tr)))) // ' s (back transform)'

  contains

    function fmt1(x) result(s)
      real, intent(in) :: x
      character(len=24) :: s
      write(s,'(f10.2)') x
    end function fmt1

    function blockof(x, y, z) result(b)
      real(dp), intent(in) :: x, y, z
      integer :: b, bx, by, bz
      bx = min(nsbx, max(1, int((x - xmin + 0.5d0*dx)/bsx) + 1))
      by = min(nsby, max(1, int((y - ymin + 0.5d0*dy)/bsy) + 1))
      bz = min(nsbz, max(1, int((z - zmin + 0.5d0*dz)/bsz) + 1))
      b = bx + nsbx*(by - 1 + nsby*(bz - 1))
    end function blockof

    subroutine krige_node(i, j, k, nsum, nfail, nnone, nexact)
      integer, intent(in) :: i, j, k
      integer(8), intent(inout) :: nsum
      integer, intent(inout) :: nfail, nnone, nexact
      real(dp) :: x, y, z, hx, hy, hz, d2, rx, ry, rz
      real(dp) :: cd(ndmax), a(ndmax, ndmax), bvec(ndmax), wgt(ndmax)
      integer :: cid(ndmax), nn, bx, by, bz, ib, jb, kb, blk, p, q, m, l, jj
      real(dp) :: s, zk, sks, ridge
      logical :: ok

      x = xmin + (i-1)*dx; y = ymin + (j-1)*dy; z = zmin + (k-1)*dz
      bx = min(nsbx, max(1, int((x - xmin + 0.5d0*dx)/bsx) + 1))
      by = min(nsby, max(1, int((y - ymin + 0.5d0*dy)/bsy) + 1))
      bz = min(nsbz, max(1, int((z - zmin + 0.5d0*dz)/bsz) + 1))
      nn = 0
      do kb = max(1, bz - mbz), min(nsbz, bz + mbz)
        do jb = max(1, by - mby), min(nsby, by + mby)
          do ib = max(1, bx - mbx), min(nsbx, bx + mbx)
            blk = ib + nsbx*(jb - 1 + nsby*(kb - 1))
            do p = sbstart(blk), sbstart(blk+1) - 1
              q = sbidx(p)
              hx = dxs(q) - x; hy = dys(q) - y; hz = dzs(q) - z
              rx = srch%rot(1,1)*hx + srch%rot(1,2)*hy + srch%rot(1,3)*hz
              ry = srch%rot(2,1)*hx + srch%rot(2,2)*hy + srch%rot(2,3)*hz
              rz = srch%rot(3,1)*hx + srch%rot(3,2)*hy + srch%rot(3,3)*hz
              d2 = (rx*rx + ry*ry + rz*rz)/(rmax*rmax)
              if (d2 > 1.0d0) cycle
              if (nn < ndmax) then
                nn = nn + 1; m = nn
              else if (d2 < cd(nn)) then
                m = nn
              else
                cycle
              end if
              do while (m > 1)
                if (cd(m-1) <= d2) exit
                cd(m) = cd(m-1); cid(m) = cid(m-1); m = m - 1
              end do
              cd(m) = d2; cid(m) = q
            end do
          end do
        end do
      end do
      nsum = nsum + nn
      if (nn == 0) then
        nnone = nnone + 1
        return
      end if
      if (cd(1) == 0.0d0) then                       ! the node itself is a data node
        if (ix(cid(1)) == i .and. iy(cid(1)) == j .and. iz(cid(1)) == k) then
          do l = 1, nreal
            res(i,j,k,l) = dg(cid(1))
          end do
          nexact = nexact + 1
          return
        end if
      end if
      ! ---- simple kriging system (covariance form) ----------------------
      do p = 1, nn
        do q = 1, p
          a(p,q) = covariance(st, nst, nugget, dxs(cid(p)) - dxs(cid(q)), &
                              dys(cid(p)) - dys(cid(q)), dzs(cid(p)) - dzs(cid(q)))
          a(q,p) = a(p,q)
        end do
        bvec(p) = covariance(st, nst, nugget, dxs(cid(p)) - x, dys(cid(p)) - y, dzs(cid(p)) - z)
      end do
      ridge = 1.0d-10*c0
      ok = .false.
      do jj = 1, 3
        call chol_solve(a, bvec, wgt, nn, ridge, ok)
        if (ok) exit
        ridge = ridge*1.0d3
      end do
      if (.not. ok) then
        nfail = nfail + 1
        return
      end if
      if (jj > 1) nfail = nfail + 1
      zk = 0.0d0
      do p = 1, nn
        zk = zk + wgt(p)*dg(cid(p))
      end do
      do l = 1, nreal
        sks = 0.0d0
        do p = 1, nn
          sks = sks + wgt(p)*simd(cid(p), l)
        end do
        res(i,j,k,l) = zk + res(i,j,k,l) - sks
      end do
    end subroutine krige_node

    subroutine chol_solve(a, bvec, wgt, n, ridge, ok)     ! solves A w = b, A SPD (copy of a is used)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(ndmax, ndmax), bvec(ndmax), ridge
      real(dp), intent(out) :: wgt(ndmax)
      logical, intent(out) :: ok
      real(dp) :: lmat(ndmax, ndmax), s
      integer :: p, q, m
      ok = .true.
      do p = 1, n
        do q = 1, p
          s = a(p,q)
          if (p == q) s = s + ridge
          do m = 1, q-1
            s = s - lmat(p,m)*lmat(q,m)
          end do
          if (p == q) then
            if (s <= 0.0d0) then
              ok = .false.; return
            end if
            lmat(p,p) = sqrt(s)
          else
            lmat(p,q) = s/lmat(q,q)
          end if
        end do
      end do
      do p = 1, n                         ! forward substitution L y = b
        s = bvec(p)
        do m = 1, p-1
          s = s - lmat(p,m)*wgt(m)
        end do
        wgt(p) = s/lmat(p,p)
      end do
      do p = n, 1, -1                     ! back substitution L^T w = y
        s = wgt(p)
        do m = p+1, n
          s = s - lmat(m,p)*wgt(m)
        end do
        wgt(p) = s/lmat(p,p)
      end do
    end subroutine chol_solve

  end subroutine condition_field

end module cond_mod

program specsim3d
  use sim3d_mod
  use cond_mod
  implicit none
  character(len=256) :: parfile, outbase, fname
  character(len=256) :: hline
  logical :: hdrfound
  integer :: ihdr, ic
  integer :: nst, nx, ny, nz, iseed, nreal, ipad
  integer :: npx, npy, npz, ntot, s, i, j, k, r, ip, jp, kp, lin, linp
  integer :: rx, ry, rz, t0, t1, tr
  real(dp) :: xmin, ymin, zmin, dx, dy, dz, nugget
  real(dp) :: ex, ey, ez, u, phi, a, c0, sneg, spos, maximag, maxreal
  real(dp) :: zmean, zvar
  type(structure) :: st(maxst)
  complex(dp), allocatable :: c(:,:,:)
  real(dp), allocatable :: amp(:,:,:), res(:,:,:,:)
  logical :: fexist, binout, vtkall, dovario, docond
  logical :: legacy, outgs, outvtk1a, outvtk1b, outvti
  integer :: nout, icond, ios, io, icode
  character(len=256) :: varfile, condfile
  integer :: narg, iarg, tsv
  real(4), allocatable :: sp(:,:,:)
  character(len=64) :: arg
  integer :: tsp, tsm, ncli

  call system_clock(t0, tr)
  parfile = 'SPECSIM3D.PAR'
  if (command_argument_count() >= 1) call get_command_argument(1, parfile)
  binout = .false.; vtkall = .false.; dovario = .false.; varfile = ''; docond = .false.; condfile = ''; ncli = 0
  narg = command_argument_count()
  iarg = 2
  do while (iarg <= narg)
    call get_command_argument(iarg, arg)
    select case (trim(arg))
    case ('--bin')
      binout = .true.; ncli = ncli + 1
    case ('--vtk')
      vtkall = .true.; ncli = ncli + 1
    case ('--vario')
      iarg = iarg + 1
      if (iarg > narg) stop '--vario needs a variogram parameter file'
      call get_command_argument(iarg, varfile)
      dovario = .true.
    case ('--cond')
      iarg = iarg + 1
      if (iarg > narg) stop '--cond needs a conditioning parameter file'
      call get_command_argument(iarg, condfile)
      docond = .true.; ncli = ncli + 1
    case default
      write(*,*) 'unknown option: ', trim(arg); stop 1
    end select
    iarg = iarg + 1
  end do
  inquire(file=trim(parfile), exist=fexist)
  if (.not. fexist) then
    write(*,*) 'Parameter file not found: ', trim(parfile); stop 1
  end if
  open(10, file=trim(parfile), status='old', action='read')
  ! Optional GSLIB-style header: any text up to a line containing 'START OF PARAMETERS:'.
  ! Files without that line (older format) are read from the first line.
  hdrfound = .false.
  do ihdr = 1, 30
    read(10,'(a)',iostat=ios) hline
    if (ios /= 0) exit
    do ic = 1, len_trim(hline)
      if (hline(ic:ic) >= 'a' .and. hline(ic:ic) <= 'z') hline(ic:ic) = achar(iachar(hline(ic:ic)) - 32)
    end do
    if (index(hline, 'START OF PARAMETERS') > 0) then
      hdrfound = .true.; exit
    end if
  end do
  if (.not. hdrfound) rewind(10)
  read(10,*) xmin, ymin, zmin
  read(10,*) dx, dy, dz
  read(10,*) nx, ny, nz
  read(10,*) nugget, nst
  if (nst < 1 .or. nst > maxst) stop 'nstruct out of range'
  do s = 1, nst
    read(10,*) st(s)%itype, st(s)%sill, st(s)%ahmax, st(s)%ahmin, &
               st(s)%avert, st(s)%az, st(s)%dip, st(s)%plg
    if (st(s)%ahmax <= 0 .or. st(s)%ahmin <= 0 .or. st(s)%avert <= 0) &
      stop 'ranges must be positive'
    call setrot(st(s))
  end do
  read(10,*) iseed, nreal, ipad
  read(10,'(a)') outbase
  outbase = adjustl(outbase)
  ! ---- optional extension of the parameter file: outputs and conditioning ----
  ! Files that end after the base name use the legacy convention (--bin, --vtk, --cond).
  outgs = .false.; outvtk1a = .false.; outvtk1b = .false.; outvti = .false.
  nout = 0; legacy = .false.
  read(10,*,iostat=ios) nout
  if (ios < 0) then
    legacy = .true.
    close(10)
  else if (ios > 0) then
    write(*,*) 'parameter file: cannot read the number of outputs (nout)'; stop 1
  end if
  if (.not. legacy) then
    binout = .false.; vtkall = .false.; docond = .false.      ! the parameter file decides
    if (nout < 1 .or. nout > 9) then
      write(*,*) 'parameter file: nout must be in 1..9'; stop 1
    end if
    do io = 1, nout
      read(10,*,iostat=ios) icode
      if (ios /= 0) then
        write(*,*) 'parameter file: missing output code'; stop 1
      end if
      select case (icode)
      case (1); outgs = .true.
      case (2); vtkall = .true.
      case (3); binout = .true.
      case (4); outvti = .true.
      case default
        write(*,*) 'parameter file: output code must be 1 (GSLIB), 2 (VTK), 3 (binary) or 4 (VTI)'; stop 1
      end select
    end do
    read(10,*,iostat=ios) icond
    if (ios /= 0) then
      write(*,*) 'parameter file: missing conditioning flag (icond)'; stop 1
    end if
    if (icond < 0 .or. icond > 1) then
      write(*,*) 'parameter file: icond must be 0 or 1'; stop 1
    end if
    docond = (icond == 1)
    if (ncli > 0) write(*,*) 'Note: --bin, --vtk and --cond are ignored: outputs and conditioning are set in the parameter file'
  else
    if (vtkall) then
      continue
    else if (binout) then
      outvtk1b = .true.
    else
      outgs = .true.; outvtk1a = .true.
    end if
  end if

  ! ---- periodic grid size ----------------------------------------
  npx = nextpow2(nx); npy = nextpow2(ny); npz = nextpow2(nz)
  if (ipad == 1) then
    ex = 0; ey = 0; ez = 0
    do s = 1, nst
      ex = max(ex, cutoff(st(s)%itype)*sqrt(sum(st(s)%axes(:,1)**2))/dx)
      ey = max(ey, cutoff(st(s)%itype)*sqrt(sum(st(s)%axes(:,2)**2))/dy)
      ez = max(ez, cutoff(st(s)%itype)*sqrt(sum(st(s)%axes(:,3)**2))/dz)
    end do
    rx = ceiling(ex); ry = ceiling(ey); rz = ceiling(ez)
    npx = nextpow2(max(nx + rx, 2*rx))
    npy = nextpow2(max(ny + ry, 2*ry))
    if (nz > 1) then
      npz = nextpow2(max(nz + rz, 2*rz))
    else
      npz = 1
    end if
  end if
  ntot = npx*npy*npz
  if (real(ntot, dp) > 2.0d0**27) stop 'periodic grid too large (>2^27 nodes)'
  write(*,'(a,3i6,a,3i6)') ' SPECSIM3D  grid:', nx, ny, nz, '   periodic grid:', npx, npy, npz

  ! ---- covariance on the periodic grid ---------------------------
  allocate(c(npx,npy,npz), amp(npx,npy,npz))
!$omp parallel do collapse(2) private(i, j, k, ip, jp, kp) schedule(static)
  do k = 1, npz
    do j = 1, npy
      do i = 1, npx
        kp = k-1; if (kp > npz/2) kp = kp - npz
        jp = j-1; if (jp > npy/2) jp = jp - npy
        ip = i-1; if (ip > npx/2) ip = ip - npx
        c(i,j,k) = cmplx(covariance(st, nst, nugget, ip*dx, jp*dy, kp*dz), 0.0d0, dp)
      end do
    end do
  end do
!$omp end parallel do
  c0 = real(c(1,1,1), dp)

  ! ---- spectral density ------------------------------------------
  call fft3d(c, npx, npy, npz, -1)
  c = c/real(ntot, dp)
  maximag = 0.0d0; maxreal = 0.0d0
  sneg = 0.0d0; spos = 0.0d0
  do k = 1, npz; do j = 1, npy; do i = 1, npx
    a = real(c(i,j,k), dp)
    maximag = max(maximag, abs(aimag(c(i,j,k))))
    maxreal = max(maxreal, abs(a))
    if (a < 0.0d0) then
      sneg = sneg + a
      amp(i,j,k) = 0.0d0
    else
      spos = spos + a
      amp(i,j,k) = sqrt(a)
    end if
  end do; end do; end do
  write(*,'(a,es10.3)') ' max |Im S| / max |Re S|        : ', maximag/maxreal
  write(*,'(a,f10.5)')  ' target variance C(0)           : ', c0
  write(*,'(a,f10.5)')  ' sum of positive spectrum       : ', spos
  write(*,'(a,es10.3)') ' negative spectral mass / total : ', abs(sneg)/(spos - sneg)

  call system_clock(tsp)
  write(*,'(a,f8.2,a)') ' time covariance+spectrum: ', real(tsp-t0)/real(tr), ' s'
  ! ---- realizations ----------------------------------------------
  rng_state = int(mod(abs(iseed), 2147483646), 8) + 1_8
  do i = 1, 10; u = urand(); end do
  allocate(res(nx,ny,nz,nreal))
  do r = 1, nreal
    do k = 0, npz-1
      do j = 0, npy-1
        do i = 0, npx-1
          ip = mod(npx-i, npx); jp = mod(npy-j, npy); kp = mod(npz-k, npz)
          lin  = i  + npx*(j  + npy*k)
          linp = ip + npx*(jp + npy*kp)
          if (lin < linp) then
            phi = 2.0d0*pi*urand()
            c(i+1,j+1,k+1)    = amp(i+1,j+1,k+1)*cmplx(cos(phi),  sin(phi), dp)
            c(ip+1,jp+1,kp+1) = amp(ip+1,jp+1,kp+1)*cmplx(cos(phi), -sin(phi), dp)
          else if (lin == linp) then
            if (urand() < 0.5d0) then
              c(i+1,j+1,k+1) = cmplx( amp(i+1,j+1,k+1), 0.0d0, dp)
            else
              c(i+1,j+1,k+1) = cmplx(-amp(i+1,j+1,k+1), 0.0d0, dp)
            end if
          end if
        end do
      end do
    end do
    call fft3d(c, npx, npy, npz, +1)
    maximag = 0.0d0
    do k = 1, nz; do j = 1, ny; do i = 1, nx
      maximag = max(maximag, abs(aimag(c(i,j,k))))
      res(i,j,k,r) = real(c(i,j,k), dp)
    end do; end do; end do
    zmean = sum(res(:,:,:,r))/real(nx*ny*nz, dp)
    zvar  = sum((res(:,:,:,r)-zmean)**2)/real(nx*ny*nz, dp)
    write(*,'(a,i3,a,f9.4,a,f9.4,a,es9.2)') ' realization', r, '  mean=', zmean, &
          '  variance=', zvar, '  max|Im|=', maximag
  end do

  if (docond) call condition_field(res, nx, ny, nz, nreal, xmin, ymin, zmin, dx, dy, dz, &
                                   st, nst, nugget, trim(condfile), trim(outbase), merge(10, 0, .not. legacy))
  if (.not. legacy) close(10)
  call system_clock(tsm)
  write(*,'(a,f8.2,a)') ' time realizations         : ', real(tsm-tsp)/real(tr), ' s'
  if (dovario) then
    call vario_report(trim(varfile), trim(outbase), res, nx, ny, nz, nreal, dx, dy, dz, st, nst, nugget)
    call system_clock(tsv)
    write(*,'(a,f8.2,a)') ' time variograms           : ', real(tsv-tsm)/real(tr), ' s'
  end if
  ! ---- output -----------------------------------------------------
  ! outgs  : GSLIB grid <base>.out        outvtk1a : ASCII VTK of realization 1 (legacy default)
  ! binout : raw float32 <base>.bin       outvtk1b : binary VTK of realization 1 (legacy --bin)
  ! vtkall : binary VTK <base>_<r>.vtk for every realization
  if (binout .or. vtkall .or. outvtk1b) then
    allocate(sp(nx,ny,nz))
    if (binout) then            ! raw float32, realizations one after the other
      fname = trim(outbase)//'.bin'
      open(11, file=trim(fname), status='replace', action='write', access='stream', form='unformatted')
    end if
    do r = 1, nreal
      do k = 1, nz; do j = 1, ny; do i = 1, nx
        sp(i,j,k) = real(res(i,j,k,r), 4)
      end do; end do; end do
      if (binout) write(11) sp
      if (vtkall .or. (outvtk1b .and. r == 1)) then
        write(fname,'(a,a,i0,a)') trim(outbase), '_', r, '.vtk'
        call write_vtk_bin(fname, sp, nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, 'SPECSIM3D realization')
      end if
    end do
    if (binout) close(11)
    deallocate(sp)
  end if
  if (outvti) then            ! one XML ImageData file with all the realizations
    fname = trim(outbase)//'.vti'
    call write_vti(fname, res, nx, ny, nz, nreal, xmin, ymin, zmin, dx, dy, dz)
  end if
  if (outgs) then
    fname = trim(outbase)//'.out'
    open(11, file=trim(fname), status='replace', action='write')
    write(11,'(a)') 'SPECSIM3D realizations'
    write(11,'(i4)') nreal
    do r = 1, nreal
      write(11,'(a,i0)') 'sim', r
    end do
    do k = 1, nz; do j = 1, ny; do i = 1, nx
      write(11,'(*(es15.7))') (res(i,j,k,r), r = 1, nreal)
    end do; end do; end do
    close(11)
  end if
  if (outvtk1a) then
    fname = trim(outbase)//'_1.vtk'
    open(12, file=trim(fname), status='replace', action='write')
    write(12,'(a)') '# vtk DataFile Version 3.0'
    write(12,'(a)') 'SPECSIM3D realization 1'
    write(12,'(a)') 'ASCII'
    write(12,'(a)') 'DATASET STRUCTURED_POINTS'
    write(12,'(a,3i6)') 'DIMENSIONS', nx, ny, nz
    write(12,'(a,3es15.7)') 'ORIGIN', xmin, ymin, zmin
    write(12,'(a,3es15.7)') 'SPACING', dx, dy, dz
    write(12,'(a,i0)') 'POINT_DATA ', nx*ny*nz
    write(12,'(a)') 'SCALARS sim float 1'
    write(12,'(a)') 'LOOKUP_TABLE default'
    do k = 1, nz; do j = 1, ny
      write(12,'(*(es13.5))') (res(i,j,k,1), i = 1, nx)
    end do; end do
    close(12)
  end if

  call system_clock(t1)
  write(*,'(a,f8.2,a)') ' elapsed time: ', real(t1-t0)/real(tr), ' s'
end program specsim3d
