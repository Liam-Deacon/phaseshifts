program hartfock_main
      implicit none
      integer :: argc, iarg
      character(len=255) :: argv, input_file, log_file, output_dir
      interface
        subroutine hartfock(inpfile, logfile, outdir)
          character(len=255) :: inpfile, logfile, outdir
        end subroutine
      end interface
    
      !Set default arguments
      input_file = 'inputA'
      log_file = 'ilogA'  !Note: Log file will use stdout if zero length string given to hartfock
      output_dir = '.'    !chgden* files will be outputted to current working directory
    
      !Check if any arguments are found
      argc = command_argument_count()
      !requires compiler adhering to the Fortran 2003 standard
      !Loop over the arguments
      iarg=1
      !loop across options
      do while(iarg <= argc)
        call get_command_argument(iarg, argv) 
        !requires compiler adhering to the Fortran 2003 standard
        select case(adjustl(argv))
        case("--help", "-h")
          write(*,*) "hartfock"
          write(*,*) "--input, -i <inputX> : path to input file"
          write(*,*) "--output-dir, -o <output_path> : path to output"
          write(*,*) "Please contact Eric Shirley <eric.shirley@nist.gov>"
          write(*,*) "for queries and comments"
          stop
        case("--input", "-i")
          if (iarg+1 .le. argc) then
            call get_command_argument(iarg+1, argv)
            input_file = adjustl(argv)
            iarg = iarg + 1
          end if
        case("--output-dir", "-o")
          if (iarg+1 .le. argc) then
            call get_command_argument(iarg+1, argv)
            output_dir = adjustl(argv)
            iarg = iarg + 1
          end if
        case("--log", "-l")
          if (iarg+1 .le. argc) then
            call get_command_argument(iarg+1, argv)
            log_file = adjustl(argv)
            iarg = iarg + 1
          end if
        end select
        iarg = iarg + 1
      end do

      call hartfock(input_file,log_file,output_dir) ! log_file is optional, and defaults to "ilogA"
      stop
end program hartfock_main
