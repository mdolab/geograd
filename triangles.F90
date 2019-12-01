module triangles
    implicit none
    real(kind=8), private :: dummyreal
    integer, parameter :: realType = kind(dummyreal)
    real(kind=realType) :: one = 1.0
    real(kind=realType) :: zero = 0.0
    contains

    subroutine dot_prod(d, v, w)
        implicit none
        real(kind=realType), dimension(3), intent(in) :: v, w
        real(kind=realType), intent(out) :: d
        d = v(1)*w(1)+v(2)*w(2)+v(3)*w(3)
    end subroutine dot_prod


    subroutine point_tri(a, b, c, p, dsquared)
        implicit none
        real(kind=realType), dimension(3), intent(in) :: a, b, c, p
        real(kind=realType), intent(out) :: dsquared

        real(kind=realType), dimension(3) :: ab, ac, ap, bp, cp, closepoint, diff, dummydiff
        real(kind=realType) :: d1, d2, d3, d4, d5, d6
        real(kind=realType) :: v, vc, vb, va, w, denom

        ab = b - a
        ac = c - a
        ap = p - a
        call dot_prod(d1, ab, ap)
        call dot_prod(d2, ac, ap)
        if (d1 <= zero .AND. d2 <= zero) then
            closepoint = a  ! barycentric 1, 0, 0
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return 
        end if

        ! check if P in vertex region outside B
        bp = p - b
        call dot_prod(d3, ab, bp)
        call dot_prod(d4, ac, bp)
        if (d3 >= zero .AND. d4 <= d3) then
            closepoint = b ! barycentric 0, 1, 0
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return
        end if

        ! check if P in edge region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*d2
        if (vc <= zero .AND. d1 >= zero .AND. d3 <= zero) then
            v = d1 / (d1 - d3)
            closepoint = a + v * ab ! barycentric coordinates (1-v,v,0)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return 
        end if
        
        ! Check if P in vertex region C
        cp = p - c
        call dot_prod(d5, ab, cp)
        call dot_prod(d6, ac, cp)
        if (d6 >= zero .AND. d5 <= d6) then
            closepoint = c ! barycentric coordinates (0,0,1)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return
        end if

        ! check if P in edge region of AC, if so, return proj(P,AC)
        vb = d5*d2 - d1*d6
        if (vb <= zero .AND. d2 >= zero .AND. d6 <=zero) then
            w = d2 / (d2 - d6)
            closepoint = a + w * ac ! barycentric (1-w, 0, w)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return
        end if

        ! Check if P in edge region of BC, if so, return proj(P,BC)
        va = d3*d6 - d5*d4
        if (va <= zero .AND. (d4-d3) >= zero .AND. (d5-d6) >= zero) then
            w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
            closepoint = b + w * (c - b) ! barycentric (0, 1-w, w)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(dsquared, diff, dummydiff)
            return
        end if

        ! P inside face region. Compute Q through barycentric (u, v, w)
        denom = one / (va + vb + vc)
        v = vb * denom
        w = vc * denom
        closepoint = a + ab * v + ac * w
        diff = closepoint - p
        dummydiff = diff
        call dot_prod(dsquared, diff, dummydiff)
        return
        
    end subroutine point_tri

    subroutine clamp(n, min, max)
        implicit none
        real(kind=realType), intent(in) :: min, max
        real(kind=realType), intent(inout) :: n

        if (n < min) then
            n = min
        end if
        if (n > max) then
            n = max
        end if
    end subroutine clamp

    subroutine line_line(p1, q1, p2, q2, dsquared)
        implicit none
        real(kind=realType), dimension(3), intent(in) :: p1, q1, p2, q2
        real(kind=realType), intent(out) :: dsquared
        real(kind=realType), dimension(3) :: d1, d2, r, diff, c1, c2, dummy3
        real(kind=realType), parameter :: EPS = 1e-12
        real(kind=realType) :: a, b, c, e, f, s, t, denom

        d1 = q1 - p1
        d2 = q2 - p2
        r = p1 - p2
        dummy3 = d1
        call dot_prod(a, d1, dummy3)
        dummy3 = d2
        call dot_prod(e, d2, dummy3)
        call dot_prod(f, d2, r)

        if (a <= EPS .AND. e <= EPS) then
            ! both segments degenrate into points
            diff = q1 - p1
            dummy3 = diff
            call dot_prod(dsquared, diff, dummy3)
            return
        end if
        if (a <= EPS) then
            s = zero
            t = f / e
            call clamp(t, zero, one)
        else
            call dot_prod(c, d1, r)
            if (e <= EPS) then
                t = zero
                s = -c / a
                call clamp(s, zero, one)
            else
                ! General non-degenerate case
                call dot_prod(b, d1, d2)
                denom = a*e - b*b
                if (denom /= zero) then
                    s = (b*f - c*e) / denom
                    call clamp(s, zero, one)
                else
                    s = zero
                end if

                t = (b*s + f)/e
                if (t < zero) then
                    t = zero
                    s = -c / a
                    call clamp(s, zero, one)
                else if (t > one) then
                    t = one
                    s = (b - c)/a
                    call clamp(s, zero, one)
                end if
            end if
        end if
        c1 = p1 + d1 * s
        c2 = p2 + d2 * t
        diff = c2 - c1
        dummy3 = diff
        call dot_prod(dsquared, diff, dummy3)
        return

    end subroutine line_line

end module triangles