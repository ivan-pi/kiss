module kiss_bmp

   use, intrinsic :: iso_fortran_env, only: int8, int16, int32

   implicit none
   private

   public :: kiss_print_bmp

   type :: bmp_header_
      sequence
      character(len=2) :: header = 'BM'
      integer(int32) :: filesize         ! unsigned
      integer(int16) :: reserved1 = 0
      integer(int16) :: reserved2 = 0
      integer(int32) :: offset           ! unsigned
   end type

   interface bmp_header_
      module procedure bmp_header_simple
   end interface

   type :: bmp_info_header_
      sequence
      integer(int32) :: headersize = 40
      integer(int32) :: width, height
      integer(int16) :: ncolorplanes = 1
      integer(int16) :: nbitsperpixel = 24    ! RGB24
      integer(int32) :: pixelcompression = 0
      integer(int32) :: datasize
      integer(int32) :: hpr = 2835
      integer(int32) :: vpr = 2835
      integer(int32) :: npalettecolors = 0
      integer(int32) :: color_importance = 0
   end type

contains

   function bmp_header_simple(filesize,offset) result(header)
      integer(int32), intent(in) :: filesize, offset
      type(bmp_header_) :: header
      header%filesize = filesize
      header%offset = offset
   end function

   !> Write integer as byte string in little endian encoding
   pure function to_bytes_i4(val) result(str)
      !> Integer value to convert to bytes
      integer, intent(in) :: val
      !> String of bytes
      character(len=4) :: str

      str = achar(mod(val, 256**1)) // &
         & achar(mod(val, 256**2) / 256**1) // &
         & achar(mod(val, 256**3) / 256**2) // &
         & achar(val / 256**3)
   end function to_bytes_i4

   !> Write integer as byte string in little endian encoding, 2-byte truncated version
   pure function to_bytes_i2(val) result(str)
      implicit none
      !> Integer value to convert to bytes
      integer, intent(in) :: val
      !> String of bytes
      character(len=2) :: str

      str = achar(mod(val, 2**8)) // &
         & achar(mod(val, 2**16) / 2**8)
   end function to_bytes_i2

   pure function bmp_header(filesize,offset) result(header)
      integer, intent(in) :: filesize, offset
      character(len=14) :: header

      !> Header field used to identify BMP file
      header(1:2) = 'BM'
      !> The size of the BMP file in bytes
      header(3:6) = to_bytes_i4(filesize)
      !> Reserved, can be 0 if created manually
      header(7:8) = achar(0)
      !> Reserved, can be 0 if created manually
      header(9:10) = achar(0)
      !> The offset, i.e. starting address, of the byte where the bitmap image
      !> data (pixel array) can be found
      header(11:14) = to_bytes_i4(offset)

   end function

   subroutine apply_colormap(value,r,g,b)
      real, intent(in) :: value
      real, intent(out) :: r, g, b

      ! For good visibility use a color scheme
      if (value <= 1./8.) then
         r = 0.
         g = 0.
         b = interp(value,-1./8.,1./8.,0.,1.)
      else if (value <= 3./8.) then
         r = 0.
         g = interp(value,1./8.,3./8.,0.,1.)
         b = 1.
      else if (value <= 5./8.) then
         r = interp(value,3./8.,5./8.,0.,1.)
         g = 1.
         b = interp(value,3./8.,5./8.,1.,0.)
      else if (value <= 7./8.) then
         r = 1.
         g = interp(value,5./8.,7./8.,1.,0.)
         b = 0.
      else
         r = interp(value,7./8.,9./8.,1.,0.)
         g = 0.
         b = 0.
      end if

   contains

      pure function interp(x,x1,x2,y1,y2) result(y)
         real, intent(in) :: x, x1, x2, y1, y2
         real :: y
         y = ((y2 - y1)*x + x2*y1 - x1*y2) / (x2 - x1)
      end function

   end subroutine

   subroutine kiss_print_bmp(data,minv,maxv,filename)

      real, intent(in) :: data(:,:)
      real, intent(in) :: minv, maxv
      character(len=*), intent(in) :: filename

      character(len=40) :: bmp_infoheader ! (DIB Header)

      integer :: filesize, datasize, w, h, padding
      integer :: unit, stat, ix, iy, pos
      integer(int8), allocatable :: img(:)
      real :: value, r, g, b

      w = size(data,2)
      h = size(data,1)

      padding = mod(4 - mod(w*3,4),4)
      !print *, "padding = ", padding

      datasize = (3*w + padding)*h
      !print *, "datasize = ", datasize

      filesize = 54 + datasize

      ! --- BITMAPINFOHEADER ---

      bmp_infoheader(1:4) = to_bytes_i4(40)  ! Header size in bytes
      bmp_infoheader(5:8) = to_bytes_i4(w)   ! Width
      bmp_infoheader(9:12) = to_bytes_i4(h)  ! Height
      bmp_infoheader(13:14) = to_bytes_i2(1) ! Number of color planes

      bmp_infoheader(15:16) = to_bytes_i2(24) ! Number of bits per pixel

      bmp_infoheader(17:20) = to_bytes_i4(0)  ! Pixel array compresion, 0 = none

      !> Size of the raw bitmap data (including padding)
      bmp_infoheader(21:24) = to_bytes_i4(datasize)

      bmp_infoheader(25:28) = to_bytes_i4(2835) ! Print resolution (horizontal)
      bmp_infoheader(29:32) = to_bytes_i4(2835) ! Print resolution (vertical)

      bmp_infoheader(33:36) = to_bytes_i4(0) ! Number of colors in the palette
      bmp_infoheader(37:40) = to_bytes_i4(0) ! Color importance



      allocate(img(datasize))

      do ix = 1, w
         do iy = 1, h

            ! Restrain values to range [0,1]
            value = (data(iy,ix) - minv) / (maxv - minv)

            call apply_colormap(value,r,g,b)

            ! Map to interval [0,255]
            r = min(255.*r,255.)
            g = min(255.*g,255.)
            b = min(255.*b,255.)

            ! Convert to 8-bit integers and pack into image buffer
            pos = ((ix-1) + (iy-1)*w)*3 + (iy-1)*padding + 1
            img(pos + 2) = int(r,int8)
            img(pos + 1) = int(g,int8)
            img(pos + 0) = int(b,int8)

         end do
      end do

      open(newunit=unit,file=trim(filename), &
         form='unformatted', &
         access='stream', &
         iostat=stat)

      if (stat == 0) &
         write(unit, iostat=stat) bmp_header_(filesize,offset=54)

      if (stat == 0) &
         write(unit, iostat=stat) bmp_infoheader

      if (stat == 0) &
         write(unit, iostat=stat) img

      close(unit, iostat=stat)

   end subroutine

end module
