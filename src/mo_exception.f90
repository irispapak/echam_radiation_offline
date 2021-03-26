MODULE mo_exception

  USE mo_io_units, ONLY: nerr, nlog 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn, em_error, em_param, em_debug
  PUBLIC :: open_log, close_log

  INTEGER, PARAMETER :: em_none  = 0   ! normal message
  INTEGER, PARAMETER :: em_info  = 1   ! informational message
  INTEGER, PARAMETER :: em_warn  = 2   ! warning message: number of warnings counted
  INTEGER, PARAMETER :: em_error = 3   ! error message: number of errors counted
  INTEGER, PARAMETER :: em_param = 4   ! report parameter value
  INTEGER, PARAMETER :: em_debug = 5   ! debugging message

  CHARACTER(len=256) :: message_text = ''         !++mgs

  LOGICAL :: l_debug = .FALSE.
  LOGICAL :: l_log   = .FALSE.

  INTEGER :: number_of_warnings  = 0
  INTEGER :: number_of_errors    = 0

CONTAINS

  SUBROUTINE finish (name, text, exit_no)

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text
    INTEGER,          INTENT(in), OPTIONAL :: exit_no

    INTEGER           :: iexit
    
    iexit=1

    IF (PRESENT(text)) THEN
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a,a,a)') 'FATAL ERROR in ', TRIM(name), ': ', TRIM(text)
      ELSE
        WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
      ENDIF
      IF (l_log) WRITE (nlog,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a)') 'FATAL ERROR in ', TRIM(name)
      ELSE
        WRITE (nerr,'(1x,a)') TRIM(name)
      ENDIF
      IF (l_log) WRITE (nlog,'(a,a)') TRIM(name)
    ENDIF

    WRITE (nerr,'(/,80("="),/)')
    IF (l_log) WRITE (nlog,'(/,80("="),/)')
    
    STOP 'mo_exception: finish ..'

  END SUBROUTINE finish

  SUBROUTINE message (name, text, out, level, all_print, adjust_right)

    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text
    INTEGER,           INTENT(in), OPTIONAL :: out
    INTEGER,           INTENT(in), OPTIONAL :: level
    LOGICAL,           INTENT(in), OPTIONAL :: all_print
    LOGICAL,           INTENT(in), OPTIONAL :: adjust_right

    INTEGER :: iout
    INTEGER :: ilevel
    LOGICAL :: lprint
    LOGICAL :: ladjustr     !++mgs renamed from ladjust to ladjustr

    CHARACTER(len=32) :: prefix

    CHARACTER(len=LEN(message_text)) :: write_text

    IF (PRESENT(all_print)) THEN
      lprint = all_print
    ELSE
      lprint = .FALSE.
    ENDIF

    IF (PRESENT(adjust_right)) THEN
      ladjustr = adjust_right 
    ELSE
      ladjustr = .FALSE.
    ENDIF

    IF (PRESENT(out)) THEN
      iout = out
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(level)) THEN
      ilevel = level
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)  ; prefix = '        '
    CASE (em_info)  ; prefix = 'INFO   :'
    CASE (em_warn)  ; prefix = 'WARNING:' ; number_of_warnings  = number_of_warnings+1
    CASE (em_error) ; prefix = 'ERROR  :' ; number_of_errors    = number_of_errors+1
    CASE (em_param) ; prefix = '---     '
    CASE (em_debug) ; prefix = 'DEBUG  :'
    END SELECT

    IF (.NOT. ladjustr) THEN
      message_text = TRIM(ADJUSTL(text))
    ENDIF
    IF (name /= '')  THEN
      message_text = TRIM(name) // ': ' // TRIM(message_text)
    ENDIF
    IF (ilevel > em_none) THEN
      message_text = TRIM(prefix) // ' ' // TRIM(message_text)
    ENDIF

    write_text = message_text

    WRITE(iout,'(1x,a)') TRIM(write_text)
    IF (l_log) WRITE(nlog,'(1x,a)') TRIM(write_text)
     
  END SUBROUTINE message

  SUBROUTINE open_log (logfile_name)

    CHARACTER(len=*), INTENT(in) :: logfile_name
    LOGICAL                      :: l_opened 

    INQUIRE (UNIT=nlog,OPENED=l_opened)

    IF (l_opened) THEN
      WRITE (message_text,'(a)') 'log file unit has been used already.'
      CALL message ('open_log', message_text)
      WRITE (message_text,'(a)') 'Close unit and reopen for log file.'
      CALL message ('open_log', message_text, level=em_warn)
      CLOSE (nlog)
    ENDIF

    OPEN (nlog,file=TRIM(logfile_name))
  
    l_log = .TRUE.

  END SUBROUTINE open_log

  SUBROUTINE close_log
    LOGICAL :: l_opened
   
    INQUIRE (UNIT=nlog,OPENED=l_opened)
    IF (l_opened) THEN
      CLOSE (nlog)
    ENDIF

    l_log = .FALSE.

  END SUBROUTINE close_log

END MODULE mo_exception
