module triangles
    implicit none
    contains

    subroutine point_tri(a, b, c, p, d)
        implicit none
        real, dimension(3), intent(in) :: a, b, c, p
        real, dimension(3), intent(out) :: d

        real, dimension(3) :: ab, ac, ap, bp
        real :: d1, d2, d3, d4
        real :: v, vc

        ab = b - a
        ac = c - a
        ap = p - a
        d1 = dot_product(ab, ap)
        d2 = dot_product(ac, ap)
        if (d1 <= 0.0 .AND. d2 <= 0.0) then
            d = a  ! barycentric 1, 0, 0
            return 
        end if

        ! check if P in vertex region outside B
        bp = p - b
        d3 = dot_product(ab, bp)
        d4 = dot_product(ac, bp)
        if (d3 >= 0.0 .AND. d4 <= d3) then
            d = b ! barycentric 0, 1, 0
            return
        end if

        ! check if P in edge region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*d2
        if (vc <= 0.0 .AND. d1 >= 0.0 .AND. d3 <= 0.0) then
            v = d1 / (d1 - d3)
            d = a + v * ab
            return 
        end if

    end subroutine point_tri

end module triangles