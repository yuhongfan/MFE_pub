Module input_namelists
       implicit none
       save
       integer :: nr, nth, nph, kunif
       real(8) :: time
       character*130 :: lbfile, b1file, b2file, b3file, &
                    v1file, v2file, v3file, dfile, efile, &
                    phfile,gridfile_in,gridfile_out,&
                    q1file,q2file,q3file
       namelist /input/ nr, nth, nph, kunif, lbfile, &
                b1file, b2file, b3file, v1file, v2file, v3file, &
                dfile, efile, phfile, gridfile_in, gridfile_out, &
                q1file, q2file, q3file, time
End Module input_namelists
