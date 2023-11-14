Subroutine main_input
       Use input_namelists
       implicit none

       open(unit=20,file="inparam_initvg",status="old",position="rewind")
       read(unit=20,nml=input)
       close(unit=20)

       return
End Subroutine
