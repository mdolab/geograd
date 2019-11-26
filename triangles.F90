module triangles
    implicit none
    contains

    subroutine point_tri(a, b, c, p, d)
        implicit none
        real, dimension(3), intent(in) :: a, b, c, p
        real, dimension(3), intent(out) :: d

        real, dimension(3) :: ab, ac, ap, bp, cp
        real :: d1, d2, d3, d4, d5, d6
        real :: v, vc, vb, va, w, denom

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
            d = a + v * ab ! barycentric coordinates (1-v,v,0)
            return 
        end if
        
        ! Check if P in vertex region C
        cp = p - c
        d5 = dot_product(ab, cp)
        d6 = dot_product(ac, cp)
        if (d6 >= 0.0 .AND. d5 <= d6) then
            d = c ! barycentric coordinates (0,0,1)
            return
        end if

        ! check if P in edge region of AC, if so, return proj(P,AC)
        vb = d5*d2 - d1*d6
        if (vb <= 0.0 .AND. d2 >= 0.0 .AND. d6 <=0.0) then
            w = d2 / (d2 - d6)
            d = a + w * ac ! barycentric (1-w, 0, w)
            return
        end if

        ! Check if P in edge region of BC, if so, return proj(P,BC)
        va = d3*d6 - d5*d4
        if (va <= 0.0 .AND. (d4-d3) >= 0.0 .AND. (d5-d6) >= 0.0) then
            w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
            d = b + w * (c - b) ! barycentric (0, 1-w, w)
            return
        end if

        ! P inside face region. Compute Q through barycentric (u, v, w)
        denom = 1.0 / (va + vb + vc)
        v = vb * denom
        w = vc * denom
        d = a + ab * v + ac * w
        return
        
    end subroutine point_tri

end module triangles