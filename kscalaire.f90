FUNCTION kscalaire(x,y,p)

  Use longr

  IMPLICIT NONE

  REAL(kind=long), INTENT(IN) :: x , y
  Integer , INTENT(IN) :: p
  REAL(kind=long)             :: kscalaire

  Select case (p)
  case(1)
     kscalaire= 1.D0
  case(2)
     kscalaire = exp(-(x+y))
  case(3)
     kscalaire = 1. /(x+y)
  case(4)
     kscalaire = x*y +1 
  end select


  RETURN
END FUNCTION kscalaire

