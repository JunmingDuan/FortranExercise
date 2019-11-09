Module TimeInfo
  Implicit None
  Integer time_info_spe(8)
  Character*10 time_info(3)
Contains
  Subroutine PrintTime
    Call date_and_time(time_info(1), time_info(2), time_info(3), time_info_spe)
    Print *, "Time: ", time_info(1), time_info(2)
  End Subroutine PrintTime
End Module TimeInfo