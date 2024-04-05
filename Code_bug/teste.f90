PROGRAM teste

  USE HDF5
  
  INTEGER(HID_T) :: out_id
  CHARACTER(LEN=3),PARAMETER :: v_name='CO2'
  CHARACTER(LEN=3),PARAMETER :: u_var='ppb'
  INTEGER,DIMENSION(2),PARAMETER :: v_shape=(/2,2/)
  REAL,DIMENSION(2,2) :: var_in
  INTEGER,PARAMETER :: irank=2
  INTEGER :: error

   INTEGER(HID_T) :: group_id(5)      ! Group identifier 

   INTEGER(HID_T) :: gr_id      ! Group identifier 
  
  INTEGER :: i,j

  CHARACTER(LEN=*),PARAMETER :: out_name='teste.h5'

  DO i=1,2
    DO j=1,2
      var_in(i,j)=2.34*i+j
    ENDDO
  ENDDO
  
   
  CALL h5_open_file(out_name,out_id)

 ! Create a group named "Chemical" in the file.
     !
  CALL h5gcreate_f(out_id, 'Chemical', group_id(1), error)
  PRINT *,group_id(1)
     ! Close the group.
     !
  CALL h5gclose_f(group_id(1), error)

  CALL h5gcreate_f(out_id, 'Meteorological', group_id(2), error)
  PRINT *,group_id(2)
     !
     ! Close the group.
     !
  CALL h5gclose_f(group_id(2), error)

  
  CALL h5gopen_f(group_id(1),'Chemical',gr_id,error)

  CALL h5gclose_f(group_id(1), error)
 
  
  

   
  CALL h5_close_file(out_id)

END PROGRAM teste