module triangles
    implicit none
    real(kind=8), parameter :: one = 1.0
    real(kind=8), parameter :: zero = 0.0
contains

    subroutine dot_prod(d, v, w)
        implicit none
        real(kind=8), dimension(3), intent(in) :: v, w
        real(kind=8), intent(out) :: d
        d = v(1) * w(1) + v(2) * w(2) + v(3) * w(3)
    end subroutine dot_prod

    subroutine cross_prod(x, u, v)
        implicit none
        real(kind=8), dimension(3), intent(in) :: u, v
        real(kind=8), dimension(3), intent(out) :: x
        x(1) = u(2) * v(3) - u(3) * v(2)
        x(2) = u(3) * v(1) - u(1) * v(3)
        x(3) = u(1) * v(2) - u(2) * v(1)
    end subroutine cross_prod

    subroutine point_tri(a, b, c, p, d)
        implicit none
        real(kind=8), dimension(3), intent(in) :: a, b, c, p
        real(kind=8), intent(out) :: d

        real(kind=8), dimension(3) :: ab, ac, ap, bp, cp, closepoint, diff, dummydiff
        real(kind=8) :: d1, d2, d3, d4, d5, d6
        real(kind=8) :: v, vc, vb, va, w, denom

        ab = b - a
        ac = c - a
        ap = p - a
        call dot_prod(d1, ab, ap)
        call dot_prod(d2, ac, ap)
        if (d1 <= zero .AND. d2 <= zero) then
            closepoint = a  ! barycentric 1, 0, 0
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
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
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
            return
        end if

        ! check if P in edge region of AB, if so return projection of P onto AB
        vc = d1 * d4 - d3 * d2
        if (vc <= zero .AND. d1 >= zero .AND. d3 <= zero) then
            v = d1 / (d1 - d3)
            closepoint = a + v * ab ! barycentric coordinates (1-v,v,0)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
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
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
            return
        end if

        ! check if P in edge region of AC, if so, return proj(P,AC)
        vb = d5 * d2 - d1 * d6
        if (vb <= zero .AND. d2 >= zero .AND. d6 <= zero) then
            w = d2 / (d2 - d6)
            closepoint = a + w * ac ! barycentric (1-w, 0, w)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
            return
        end if

        ! Check if P in edge region of BC, if so, return proj(P,BC)
        va = d3 * d6 - d5 * d4
        if (va <= zero .AND. (d4 - d3) >= zero .AND. (d5 - d6) >= zero) then
            w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
            closepoint = b + w * (c - b) ! barycentric (0, 1-w, w)
            diff = closepoint - p
            dummydiff = diff
            call dot_prod(d, diff, dummydiff)
            d = sqrt(d)
            return
        end if

        ! P inside face region. Compute Q through barycentric (u, v, w)
        denom = one / (va + vb + vc)
        v = vb * denom
        w = vc * denom
        closepoint = a + ab * v + ac * w
        diff = closepoint - p
        dummydiff = diff
        call dot_prod(d, diff, dummydiff)
        d = sqrt(d)
        return

    end subroutine point_tri

    subroutine clamp(n, min, max)
        implicit none
        real(kind=8), intent(in) :: min, max
        real(kind=8), intent(inout) :: n

        if (n < min) then
            n = min
        end if
        if (n > max) then
            n = max
        end if
    end subroutine clamp

    subroutine line_line(p1, q1, p2, q2, d)
        implicit none
        real(kind=8), dimension(3), intent(in) :: p1, q1, p2, q2
        real(kind=8), intent(out) :: d
        real(kind=8), dimension(3) :: d1, d2, r, diff, c1, c2, dummy3
        real(kind=8), parameter :: EPS = 1e-12
        real(kind=8) :: a, b, c, e, f, s, t, denom

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
            call dot_prod(d, diff, dummy3)
            d = sqrt(d)
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
                denom = a * e - b * b
                if (denom /= zero) then
                    s = (b * f - c * e) / denom
                    call clamp(s, zero, one)
                else
                    s = zero
                end if

                t = (b * s + f) / e
                if (t < zero) then
                    t = zero
                    s = -c / a
                    call clamp(s, zero, one)
                else if (t > one) then
                    t = one
                    s = (b - c) / a
                    call clamp(s, zero, one)
                end if
            end if
        end if
        c1 = p1 + d1 * s
        c2 = p2 + d2 * t
        diff = c2 - c1
        dummy3 = diff
        call dot_prod(d, diff, dummy3)
        d = sqrt(d)
        return

    end subroutine line_line

    subroutine maxloc3(a, maxind)
        implicit none
        real(kind=8), dimension(3), intent(in) :: a
        integer, intent(out) :: maxind

        if (a(1) > a(2)) then
            if (a(1) > a(3)) then
                maxind = 1
                return
            else
                maxind = 3
                return
            end if
        elseif (a(2) > a(3)) then
            maxind = 2
            return
        else
            maxind = 3
            return
        end if
    end subroutine maxloc3

    subroutine intersect(a1, b1, c1, a2, b2, c2, length)
        implicit none
        real(kind=8), dimension(3), intent(in) :: a1, b1, c1, a2, b2, c2
        real(kind=8), intent(out) :: length
        real(kind=8), dimension(3) :: n1, n2, d, absd
        real(kind=8), parameter :: EPS = 1e-12
        real(kind=8) :: d1, da1, db1, dc1, pa1, pb1, pc1, t11, t12, t1high, t1low, dt
        real(kind=8) :: d2, da2, db2, dc2, pa2, pb2, pc2, t21, t22, t2high, t2low
        integer :: lone_vertex_1, lone_vertex_2, maxind
        call cross_prod(n2, (b2 - a2), (c2 - a2))
        call dot_prod(d2, n2, a2)
        d2 = -d2

        call dot_prod(da1, n2, a1)
        da1 = da1 + d2
        call dot_prod(db1, n2, b1)
        db1 = db1 + d2
        call dot_prod(dc1, n2, c1)
        dc1 = dc1 + d2

        if ((da1 .GE. zero) .and. (db1 .GE. zero) .and. (dc1 .GE. zero)) then
            length = 0.0
            return
        elseif ((da1 .LE. zero) .and. (db1 .LE. zero) .and. (dc1 .LE. zero)) then
            length = 0.0
            return
        end if

        ! general case
        call cross_prod(n1, (b1 - a1), (c1 - a1))
        call dot_prod(d1, n1, a1)
        d1 = -d1

        call dot_prod(da2, n1, a2)
        da2 = da2 + d1
        call dot_prod(db2, n1, b2)
        db2 = db2 + d1
        call dot_prod(dc2, n1, c2)
        dc2 = dc2 + d1

        if ((da2 .GE. zero) .and. (db2 .GE. zero) .and. (dc2 .GE. zero)) then
            length = 0.0
            return
        elseif ((da2 .LE. zero) .and. (db2 .LE. zero) .and. (dc2 .LE. zero)) then
            length = 0.0
            return
        end if

        call cross_prod(d, n1, n2)
        ! absd = abs(d)
        ! call maxloc3(absd, maxind)
        ! pa1 = a1(maxind)
        ! pb1 = b1(maxind)
        ! pc1 = c1(maxind)
        call dot_prod(pa1, d, a1)
        call dot_prod(pb1, d, b1)
        call dot_prod(pc1, d, c1)

        ! need to figure out which vertex is by itself
        if (da1 > 0) then
            if (db1 > 0) then
                lone_vertex_1 = 3
            elseif (dc1 > 0) then
                lone_vertex_1 = 2
            else
                lone_vertex_1 = 1
            end if
        else
            if (db1 < 0) then
                lone_vertex_1 = 3
            elseif (dc1 < 0) then
                lone_vertex_1 = 2
            else
                lone_vertex_1 = 1
            end if
        end if

        if (lone_vertex_1 == 1) then
            t11 = pb1 + (pa1 - pb1) * db1 / (db1 - da1)
            t12 = pc1 + (pa1 - pc1) * dc1 / (dc1 - da1)
        elseif (lone_vertex_1 == 2) then
            t11 = pa1 + (pb1 - pa1) * da1 / (da1 - db1)
            t12 = pc1 + (pb1 - pc1) * dc1 / (dc1 - db1)
        else
            t11 = pa1 + (pc1 - pa1) * da1 / (da1 - dc1)
            t12 = pb1 + (pc1 - pb1) * db1 / (db1 - dc1)
        end if

        call dot_prod(pa2, d, a2)
        call dot_prod(pb2, d, b2)
        call dot_prod(pc2, d, c2)

        ! need to figure out which vertex is by itself
        if (da2 > 0) then
            if (db2 > 0) then
                lone_vertex_2 = 3
            elseif (dc2 > 0) then
                lone_vertex_2 = 2
            else
                lone_vertex_2 = 1
            end if
        else
            if (db2 < 0) then
                lone_vertex_2 = 3
            elseif (dc2 < 0) then
                lone_vertex_2 = 2
            else
                lone_vertex_2 = 1
            end if
        end if

        if (lone_vertex_2 == 1) then
            t21 = pb2 + (pa2 - pb2) * db2 / (db2 - da2)
            t22 = pc2 + (pa2 - pc2) * dc2 / (dc2 - da2)
        elseif (lone_vertex_2 == 2) then
            t21 = pa2 + (pb2 - pa2) * da2 / (da2 - db2)
            t22 = pc2 + (pb2 - pc2) * dc2 / (dc2 - db2)
        else
            t21 = pa2 + (pc2 - pa2) * da2 / (da2 - dc2)
            t22 = pb2 + (pc2 - pb2) * db2 / (db2 - dc2)
        end if

        if (t11 > t12) then
            t1high = t11
            t1low = t12
        else
            t1high = t12
            t1low = t11
        end if
        if (t21 > t22) then
            t2high = t21
            t2low = t22
        else
            t2high = t22
            t2low = t21
        end if
        ! check if intervals overlap
        if ((t1high < t2low) .OR. (t2high < t1low)) then
            ! no overlap
            length = 0.0
            return
        else
            dt = min(t1high, t2high) - max(t1low, t2low)
            length = dt / sqrt((d(1)**2 + d(2)**2 + d(3)**2))
            return
        end if
    end subroutine intersect

end module triangles
