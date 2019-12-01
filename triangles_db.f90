!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
!
MODULE TRIANGLES_DB
  IMPLICIT NONE
  REAL(kind=8) :: one=1.0
  REAL(kind=8) :: oned=0.0_8
  REAL(kind=8) :: oneb=0.0_8
  REAL(kind=8) :: zero=0.0
  REAL(kind=8) :: zerod=0.0_8
  REAL(kind=8) :: zerob=0.0_8

CONTAINS
!  Differentiation of point_tri in forward (tangent) mode:
!   variations   of useful results: dsquared
!   with respect to varying inputs: one p a b c
!   RW status of diff variables: one:in dsquared:out p:in a:in
!                b:in c:in
  SUBROUTINE POINT_TRI_D(a, ad, b, bd, c, cd, p, pd, dsquared, dsquaredd&
& )
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a, b, c, p
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: ad, bd, cd, pd
    REAL(kind=8), INTENT(OUT) :: dsquared
    REAL(kind=8), INTENT(OUT) :: dsquaredd
    REAL(kind=8), DIMENSION(3) :: ab, ac, ap, bp, cp, closepoint, diff, &
&   dummydiff
    REAL(kind=8), DIMENSION(3) :: abd, acd, apd, bpd, cpd, closepointd, &
&   diffd, dummydiffd
    REAL(kind=8) :: d1, d2, d3, d4, d5, d6
    REAL(kind=8) :: d1d, d2d, d3d, d4d, d5d, d6d
    REAL(kind=8) :: v, vc, vb, va, w, denom
    REAL(kind=8) :: vd, vcd, vbd, vad, wd, denomd
    abd = bd - ad
    ab = b - a
    acd = cd - ad
    ac = c - a
    apd = pd - ad
    ap = p - a
    CALL DOT_PROD_D(d1, d1d, ab, abd, ap, apd)
    CALL DOT_PROD_D(d2, d2d, ac, acd, ap, apd)
    IF (d1 .LE. zero .AND. d2 .LE. zero) THEN
! barycentric 1, 0, 0
      closepointd = ad
      closepoint = a
      diffd = closepointd - pd
      diff = closepoint - p
      dummydiffd = diffd
      dummydiff = diff
      CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummydiff, &
&               dummydiffd)
      RETURN
    ELSE
! check if P in vertex region outside B
      bpd = pd - bd
      bp = p - b
      CALL DOT_PROD_D(d3, d3d, ab, abd, bp, bpd)
      CALL DOT_PROD_D(d4, d4d, ac, acd, bp, bpd)
      IF (d3 .GE. zero .AND. d4 .LE. d3) THEN
! barycentric 0, 1, 0
        closepointd = bd
        closepoint = b
        diffd = closepointd - pd
        diff = closepoint - p
        dummydiffd = diffd
        dummydiff = diff
        CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummydiff, &
&                 dummydiffd)
        RETURN
      ELSE
! check if P in edge region of AB, if so return projection of P onto AB
        vcd = d1d*d4 + d1*d4d - d3d*d2 - d3*d2d
        vc = d1*d4 - d3*d2
        IF (vc .LE. zero .AND. d1 .GE. zero .AND. d3 .LE. zero) THEN
          vd = (d1d*(d1-d3)-d1*(d1d-d3d))/(d1-d3)**2
          v = d1/(d1-d3)
! barycentric coordinates (1-v,v,0)
          closepointd = ad + vd*ab + v*abd
          closepoint = a + v*ab
          diffd = closepointd - pd
          diff = closepoint - p
          dummydiffd = diffd
          dummydiff = diff
          CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummydiff, &
&                   dummydiffd)
          RETURN
        ELSE
! Check if P in vertex region C
          cpd = pd - cd
          cp = p - c
          CALL DOT_PROD_D(d5, d5d, ab, abd, cp, cpd)
          CALL DOT_PROD_D(d6, d6d, ac, acd, cp, cpd)
          IF (d6 .GE. zero .AND. d5 .LE. d6) THEN
! barycentric coordinates (0,0,1)
            closepointd = cd
            closepoint = c
            diffd = closepointd - pd
            diff = closepoint - p
            dummydiffd = diffd
            dummydiff = diff
            CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummydiff&
&                     , dummydiffd)
            RETURN
          ELSE
! check if P in edge region of AC, if so, return proj(P,AC)
            vbd = d5d*d2 + d5*d2d - d1d*d6 - d1*d6d
            vb = d5*d2 - d1*d6
            IF (vb .LE. zero .AND. d2 .GE. zero .AND. d6 .LE. zero) THEN
              wd = (d2d*(d2-d6)-d2*(d2d-d6d))/(d2-d6)**2
              w = d2/(d2-d6)
! barycentric (1-w, 0, w)
              closepointd = ad + wd*ac + w*acd
              closepoint = a + w*ac
              diffd = closepointd - pd
              diff = closepoint - p
              dummydiffd = diffd
              dummydiff = diff
              CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, &
&                       dummydiff, dummydiffd)
              RETURN
            ELSE
! Check if P in edge region of BC, if so, return proj(P,BC)
              vad = d3d*d6 + d3*d6d - d5d*d4 - d5*d4d
              va = d3*d6 - d5*d4
              IF (va .LE. zero .AND. d4 - d3 .GE. zero .AND. d5 - d6 &
&                 .GE. zero) THEN
                wd = ((d4d-d3d)*(d4-d3+(d5-d6))-(d4-d3)*(d4d-d3d+d5d-d6d&
&                 ))/(d4-d3+(d5-d6))**2
                w = (d4-d3)/(d4-d3+(d5-d6))
! barycentric (0, 1-w, w)
                closepointd = bd + wd*(c-b) + w*(cd-bd)
                closepoint = b + w*(c-b)
                diffd = closepointd - pd
                diff = closepoint - p
                dummydiffd = diffd
                dummydiff = diff
                CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, &
&                         dummydiff, dummydiffd)
                RETURN
              ELSE
! P inside face region. Compute Q through barycentric (u, v, w)
                denomd = (oned*(va+vb+vc)-one*(vad+vbd+vcd))/(va+vb+vc)&
&                 **2
                denom = one/(va+vb+vc)
                vd = vbd*denom + vb*denomd
                v = vb*denom
                wd = vcd*denom + vc*denomd
                w = vc*denom
                closepointd = ad + abd*v + ab*vd + acd*w + ac*wd
                closepoint = a + ab*v + ac*w
                diffd = closepointd - pd
                diff = closepoint - p
                dummydiffd = diffd
                dummydiff = diff
                CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, &
&                         dummydiff, dummydiffd)
                RETURN
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END SUBROUTINE POINT_TRI_D
!  Differentiation of point_tri in reverse (adjoint) mode:
!   gradient     of useful results: one dsquared p a b c
!   with respect to varying inputs: one dsquared p a b c
!   RW status of diff variables: one:incr dsquared:in-zero p:incr
!                a:incr b:incr c:incr
  SUBROUTINE POINT_TRI_B(a, ab0, b, bb, c, cb, p, pb, dsquared, &
&   dsquaredb)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a, b, c, p
    REAL(kind=8), DIMENSION(3) :: ab0, bb, cb, pb
    REAL(kind=8) :: dsquared
    REAL(kind=8) :: dsquaredb
    REAL(kind=8), DIMENSION(3) :: ab, ac, ap, bp, cp, closepoint, diff, &
&   dummydiff
    REAL(kind=8), DIMENSION(3) :: abb, acb, apb, bpb, cpb, closepointb, &
&   diffb, dummydiffb
    REAL(kind=8) :: d1, d2, d3, d4, d5, d6
    REAL(kind=8) :: d1b, d2b, d3b, d4b, d5b, d6b
    REAL(kind=8) :: v, vc, vb, va, w, denom
    REAL(kind=8) :: vb0, vcb, vbb, vab, wb, denomb
    REAL(kind=8) :: tempb6
    REAL(kind=8) :: tempb5
    REAL(kind=8) :: tempb4
    REAL(kind=8) :: tempb3
    REAL(kind=8) :: tempb2
    REAL(kind=8) :: tempb1
    REAL(kind=8) :: tempb0
    REAL(kind=8) :: tempb
    ab = b - a
    ac = c - a
    ap = p - a
    CALL DOT_PROD(d1, ab, ap)
    CALL DOT_PROD(d2, ac, ap)
    IF (d1 .LE. zero .AND. d2 .LE. zero) THEN
! barycentric 1, 0, 0
      closepoint = a
      diff = closepoint - p
      dummydiff = diff
      diffb = 0.0_8
      dummydiffb = 0.0_8
      CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummydiff, &
&               dummydiffb)
      diffb = diffb + dummydiffb
      closepointb = 0.0_8
      closepointb = diffb
      pb = pb - diffb
      ab0 = ab0 + closepointb
      abb = 0.0_8
      acb = 0.0_8
      d1b = 0.0_8
      d2b = 0.0_8
    ELSE
! check if P in vertex region outside B
      bp = p - b
      CALL DOT_PROD(d3, ab, bp)
      CALL DOT_PROD(d4, ac, bp)
      IF (d3 .GE. zero .AND. d4 .LE. d3) THEN
! barycentric 0, 1, 0
        closepoint = b
        diff = closepoint - p
        dummydiff = diff
        diffb = 0.0_8
        dummydiffb = 0.0_8
        CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummydiff, &
&                 dummydiffb)
        diffb = diffb + dummydiffb
        closepointb = 0.0_8
        closepointb = diffb
        pb = pb - diffb
        bb = bb + closepointb
        abb = 0.0_8
        acb = 0.0_8
        d1b = 0.0_8
        d2b = 0.0_8
        d3b = 0.0_8
        d4b = 0.0_8
      ELSE
! check if P in edge region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*d2
        IF (vc .LE. zero .AND. d1 .GE. zero .AND. d3 .LE. zero) THEN
          v = d1/(d1-d3)
! barycentric coordinates (1-v,v,0)
          closepoint = a + v*ab
          diff = closepoint - p
          dummydiff = diff
          diffb = 0.0_8
          dummydiffb = 0.0_8
          CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummydiff, &
&                   dummydiffb)
          diffb = diffb + dummydiffb
          closepointb = 0.0_8
          closepointb = diffb
          pb = pb - diffb
          abb = 0.0_8
          ab0 = ab0 + closepointb
          vb0 = SUM(ab*closepointb)
          abb = v*closepointb
          tempb = vb0/(d1-d3)
          tempb0 = -(d1*tempb/(d1-d3))
          d1b = tempb0 + tempb
          d3b = -tempb0
          acb = 0.0_8
          d2b = 0.0_8
          d4b = 0.0_8
          vcb = 0.0_8
        ELSE
! Check if P in vertex region C
          cp = p - c
          CALL DOT_PROD(d5, ab, cp)
          CALL DOT_PROD(d6, ac, cp)
          IF (d6 .GE. zero .AND. d5 .LE. d6) THEN
! barycentric coordinates (0,0,1)
            closepoint = c
            diff = closepoint - p
            dummydiff = diff
            diffb = 0.0_8
            dummydiffb = 0.0_8
            CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummydiff&
&                     , dummydiffb)
            diffb = diffb + dummydiffb
            closepointb = 0.0_8
            closepointb = diffb
            pb = pb - diffb
            cb = cb + closepointb
            abb = 0.0_8
            acb = 0.0_8
            d1b = 0.0_8
            d2b = 0.0_8
            d3b = 0.0_8
            d4b = 0.0_8
            d5b = 0.0_8
            d6b = 0.0_8
            vcb = 0.0_8
          ELSE
! check if P in edge region of AC, if so, return proj(P,AC)
            vb = d5*d2 - d1*d6
            IF (vb .LE. zero .AND. d2 .GE. zero .AND. d6 .LE. zero) THEN
              w = d2/(d2-d6)
! barycentric (1-w, 0, w)
              closepoint = a + w*ac
              diff = closepoint - p
              dummydiff = diff
              diffb = 0.0_8
              dummydiffb = 0.0_8
              CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, &
&                       dummydiff, dummydiffb)
              diffb = diffb + dummydiffb
              closepointb = 0.0_8
              closepointb = diffb
              pb = pb - diffb
              acb = 0.0_8
              ab0 = ab0 + closepointb
              wb = SUM(ac*closepointb)
              acb = w*closepointb
              tempb1 = wb/(d2-d6)
              tempb2 = -(d2*tempb1/(d2-d6))
              d2b = tempb2 + tempb1
              d6b = -tempb2
              abb = 0.0_8
              d3b = 0.0_8
              d4b = 0.0_8
              d5b = 0.0_8
              vbb = 0.0_8
              vcb = 0.0_8
            ELSE
! Check if P in edge region of BC, if so, return proj(P,BC)
              va = d3*d6 - d5*d4
              IF (va .LE. zero .AND. d4 - d3 .GE. zero .AND. d5 - d6 &
&                 .GE. zero) THEN
                w = (d4-d3)/(d4-d3+(d5-d6))
! barycentric (0, 1-w, w)
                closepoint = b + w*(c-b)
                diff = closepoint - p
                dummydiff = diff
                diffb = 0.0_8
                dummydiffb = 0.0_8
                CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, &
&                         dummydiff, dummydiffb)
                diffb = diffb + dummydiffb
                closepointb = 0.0_8
                closepointb = diffb
                pb = pb - diffb
                bb = bb + closepointb - w*closepointb
                wb = SUM((c-b)*closepointb)
                cb = cb + w*closepointb
                tempb3 = wb/(d4-d3+d5-d6)
                tempb4 = -((d4-d3)*tempb3/(d4-d3+d5-d6))
                d4b = tempb4 + tempb3
                d3b = -tempb4 - tempb3
                d5b = tempb4
                d6b = -tempb4
                abb = 0.0_8
                acb = 0.0_8
                vab = 0.0_8
                vbb = 0.0_8
                vcb = 0.0_8
              ELSE
! P inside face region. Compute Q through barycentric (u, v, w)
                denom = one/(va+vb+vc)
                v = vb*denom
                w = vc*denom
                closepoint = a + ab*v + ac*w
                diff = closepoint - p
                dummydiff = diff
                diffb = 0.0_8
                dummydiffb = 0.0_8
                CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, &
&                         dummydiff, dummydiffb)
                diffb = diffb + dummydiffb
                closepointb = 0.0_8
                closepointb = diffb
                pb = pb - diffb
                abb = 0.0_8
                acb = 0.0_8
                ab0 = ab0 + closepointb
                abb = v*closepointb
                vb0 = SUM(ab*closepointb)
                acb = w*closepointb
                wb = SUM(ac*closepointb)
                denomb = vb*vb0 + vc*wb
                tempb6 = denomb/(va+vb+vc)
                tempb5 = -(one*tempb6/(va+vb+vc))
                vcb = tempb5 + denom*wb
                vbb = tempb5 + denom*vb0
                oneb = oneb + tempb6
                vab = tempb5
                d3b = 0.0_8
                d4b = 0.0_8
                d5b = 0.0_8
                d6b = 0.0_8
              END IF
              d3b = d3b + d6*vab
              d6b = d6b + d3*vab
              d5b = d5b - d4*vab
              d4b = d4b - d5*vab
              d2b = 0.0_8
            END IF
            d5b = d5b + d2*vbb
            d2b = d2b + d5*vbb
            d1b = -(d6*vbb)
            d6b = d6b - d1*vbb
          END IF
          cpb = 0.0_8
          CALL DOT_PROD_B(d6, d6b, ac, acb, cp, cpb)
          CALL DOT_PROD_B(d5, d5b, ab, abb, cp, cpb)
          pb = pb + cpb
          cb = cb - cpb
        END IF
        d1b = d1b + d4*vcb
        d4b = d4b + d1*vcb
        d3b = d3b - d2*vcb
        d2b = d2b - d3*vcb
      END IF
      bpb = 0.0_8
      CALL DOT_PROD_B(d4, d4b, ac, acb, bp, bpb)
      CALL DOT_PROD_B(d3, d3b, ab, abb, bp, bpb)
      pb = pb + bpb
      bb = bb - bpb
    END IF
    apb = 0.0_8
    CALL DOT_PROD_B(d2, d2b, ac, acb, ap, apb)
    CALL DOT_PROD_B(d1, d1b, ab, abb, ap, apb)
    pb = pb + apb
    ab0 = ab0 - acb - abb - apb
    cb = cb + acb
    bb = bb + abb
    dsquaredb = 0.0_8
  END SUBROUTINE POINT_TRI_B
  SUBROUTINE POINT_TRI(a, b, c, p, dsquared)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a, b, c, p
    REAL(kind=8), INTENT(OUT) :: dsquared
    REAL(kind=8), DIMENSION(3) :: ab, ac, ap, bp, cp, closepoint, diff, &
&   dummydiff
    REAL(kind=8) :: d1, d2, d3, d4, d5, d6
    REAL(kind=8) :: v, vc, vb, va, w, denom
    ab = b - a
    ac = c - a
    ap = p - a
    CALL DOT_PROD(d1, ab, ap)
    CALL DOT_PROD(d2, ac, ap)
    IF (d1 .LE. zero .AND. d2 .LE. zero) THEN
! barycentric 1, 0, 0
      closepoint = a
      diff = closepoint - p
      dummydiff = diff
      CALL DOT_PROD(dsquared, diff, dummydiff)
      RETURN
    ELSE
! check if P in vertex region outside B
      bp = p - b
      CALL DOT_PROD(d3, ab, bp)
      CALL DOT_PROD(d4, ac, bp)
      IF (d3 .GE. zero .AND. d4 .LE. d3) THEN
! barycentric 0, 1, 0
        closepoint = b
        diff = closepoint - p
        dummydiff = diff
        CALL DOT_PROD(dsquared, diff, dummydiff)
        RETURN
      ELSE
! check if P in edge region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*d2
        IF (vc .LE. zero .AND. d1 .GE. zero .AND. d3 .LE. zero) THEN
          v = d1/(d1-d3)
! barycentric coordinates (1-v,v,0)
          closepoint = a + v*ab
          diff = closepoint - p
          dummydiff = diff
          CALL DOT_PROD(dsquared, diff, dummydiff)
          RETURN
        ELSE
! Check if P in vertex region C
          cp = p - c
          CALL DOT_PROD(d5, ab, cp)
          CALL DOT_PROD(d6, ac, cp)
          IF (d6 .GE. zero .AND. d5 .LE. d6) THEN
! barycentric coordinates (0,0,1)
            closepoint = c
            diff = closepoint - p
            dummydiff = diff
            CALL DOT_PROD(dsquared, diff, dummydiff)
            RETURN
          ELSE
! check if P in edge region of AC, if so, return proj(P,AC)
            vb = d5*d2 - d1*d6
            IF (vb .LE. zero .AND. d2 .GE. zero .AND. d6 .LE. zero) THEN
              w = d2/(d2-d6)
! barycentric (1-w, 0, w)
              closepoint = a + w*ac
              diff = closepoint - p
              dummydiff = diff
              CALL DOT_PROD(dsquared, diff, dummydiff)
              RETURN
            ELSE
! Check if P in edge region of BC, if so, return proj(P,BC)
              va = d3*d6 - d5*d4
              IF (va .LE. zero .AND. d4 - d3 .GE. zero .AND. d5 - d6 &
&                 .GE. zero) THEN
                w = (d4-d3)/(d4-d3+(d5-d6))
! barycentric (0, 1-w, w)
                closepoint = b + w*(c-b)
                diff = closepoint - p
                dummydiff = diff
                CALL DOT_PROD(dsquared, diff, dummydiff)
                RETURN
              ELSE
! P inside face region. Compute Q through barycentric (u, v, w)
                denom = one/(va+vb+vc)
                v = vb*denom
                w = vc*denom
                closepoint = a + ab*v + ac*w
                diff = closepoint - p
                dummydiff = diff
                CALL DOT_PROD(dsquared, diff, dummydiff)
                RETURN
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END SUBROUTINE POINT_TRI
!  Differentiation of line_line in forward (tangent) mode:
!   variations   of useful results: dsquared
!   with respect to varying inputs: zero one p1 p2 q1 q2
!   RW status of diff variables: zero:in one:in dsquared:out p1:in
!                p2:in q1:in q2:in
  SUBROUTINE LINE_LINE_D(p1, p1d, q1, q1d, p2, p2d, q2, q2d, dsquared, &
&   dsquaredd)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: p1, q1, p2, q2
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: p1d, q1d, p2d, q2d
    REAL(kind=8), INTENT(OUT) :: dsquared
    REAL(kind=8), INTENT(OUT) :: dsquaredd
    REAL(kind=8), DIMENSION(3) :: d1, d2, r, diff, c1, c2, dummy3
    REAL(kind=8), DIMENSION(3) :: d1d, d2d, rd, diffd, c1d, c2d, dummy3d
    REAL(kind=8), PARAMETER :: eps=1e-12
    REAL(kind=8) :: a, b, c, e, f, s, t, denom
    REAL(kind=8) :: ad, bd, cd, ed, fd, sd, td, denomd
    d1d = q1d - p1d
    d1 = q1 - p1
    d2d = q2d - p2d
    d2 = q2 - p2
    rd = p1d - p2d
    r = p1 - p2
    dummy3d = d1d
    dummy3 = d1
    CALL DOT_PROD_D(a, ad, d1, d1d, dummy3, dummy3d)
    dummy3d = d2d
    dummy3 = d2
    CALL DOT_PROD_D(e, ed, d2, d2d, dummy3, dummy3d)
    CALL DOT_PROD_D(f, fd, d2, d2d, r, rd)
    IF (a .LE. eps .AND. e .LE. eps) THEN
! both segments degenrate into points
      diffd = q1d - p1d
      diff = q1 - p1
      dummy3d = diffd
      dummy3 = diff
      CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummy3, dummy3d)
      RETURN
    ELSE
      IF (a .LE. eps) THEN
        sd = zerod
        s = zero
        td = (fd*e-f*ed)/e**2
        t = f/e
        CALL CLAMP_D(t, td, zero, zerod, one, oned)
      ELSE
        CALL DOT_PROD_D(c, cd, d1, d1d, r, rd)
        IF (e .LE. eps) THEN
          td = zerod
          t = zero
          sd = -((cd*a-c*ad)/a**2)
          s = -(c/a)
          CALL CLAMP_D(s, sd, zero, zerod, one, oned)
        ELSE
! General non-degenerate case
          CALL DOT_PROD_D(b, bd, d1, d1d, d2, d2d)
          denomd = ad*e + a*ed - bd*b - b*bd
          denom = a*e - b*b
          IF (denom .NE. zero) THEN
            sd = ((bd*f+b*fd-cd*e-c*ed)*denom-(b*f-c*e)*denomd)/denom**2
            s = (b*f-c*e)/denom
            CALL CLAMP_D(s, sd, zero, zerod, one, oned)
          ELSE
            sd = zerod
            s = zero
          END IF
          td = ((bd*s+b*sd+fd)*e-(b*s+f)*ed)/e**2
          t = (b*s+f)/e
          IF (t .LT. zero) THEN
            td = zerod
            t = zero
            sd = -((cd*a-c*ad)/a**2)
            s = -(c/a)
            CALL CLAMP_D(s, sd, zero, zerod, one, oned)
          ELSE IF (t .GT. one) THEN
            td = oned
            t = one
            sd = ((bd-cd)*a-(b-c)*ad)/a**2
            s = (b-c)/a
            CALL CLAMP_D(s, sd, zero, zerod, one, oned)
          END IF
        END IF
      END IF
      c1d = p1d + d1d*s + d1*sd
      c1 = p1 + d1*s
      c2d = p2d + d2d*t + d2*td
      c2 = p2 + d2*t
      diffd = c2d - c1d
      diff = c2 - c1
      dummy3d = diffd
      dummy3 = diff
      CALL DOT_PROD_D(dsquared, dsquaredd, diff, diffd, dummy3, dummy3d)
      RETURN
    END IF
  END SUBROUTINE LINE_LINE_D
!  Differentiation of line_line in reverse (adjoint) mode:
!   gradient     of useful results: zero one dsquared p1 p2 q1
!                q2
!   with respect to varying inputs: zero one dsquared p1 p2 q1
!                q2
!   RW status of diff variables: zero:incr one:incr dsquared:in-zero
!                p1:incr p2:incr q1:incr q2:incr
  SUBROUTINE LINE_LINE_B(p1, p1b, q1, q1b, p2, p2b, q2, q2b, dsquared, &
&   dsquaredb)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: p1, q1, p2, q2
    REAL(kind=8), DIMENSION(3) :: p1b, q1b, p2b, q2b
    REAL(kind=8) :: dsquared
    REAL(kind=8) :: dsquaredb
    REAL(kind=8), DIMENSION(3) :: d1, d2, r, diff, c1, c2, dummy3
    REAL(kind=8), DIMENSION(3) :: d1b, d2b, rb, diffb, c1b, c2b, dummy3b
    REAL(kind=8), PARAMETER :: eps=1e-12
    REAL(kind=8) :: a, b, c, e, f, s, t, denom
    REAL(kind=8) :: ab, bb, cb, eb, fb, sb, tb, denomb
    INTEGER :: branch
    REAL(kind=8) :: tempb1
    REAL(kind=8) :: tempb0
    REAL(kind=8) :: tempb
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    dummy3 = d1
    CALL DOT_PROD(a, d1, dummy3)
    dummy3 = d2
    CALL DOT_PROD(e, d2, dummy3)
    CALL DOT_PROD(f, d2, r)
    IF (a .LE. eps .AND. e .LE. eps) THEN
! both segments degenrate into points
      diff = q1 - p1
      dummy3 = diff
      diffb = 0.0_8
      dummy3b = 0.0_8
      CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummy3, dummy3b)
      diffb = diffb + dummy3b
      q1b = q1b + diffb
      p1b = p1b - diffb
      eb = 0.0_8
      fb = 0.0_8
      rb = 0.0_8
      d1b = 0.0_8
      d2b = 0.0_8
      ab = 0.0_8
    ELSE
      IF (a .LE. eps) THEN
        s = zero
        t = f/e
        CALL PUSHREAL8(t)
        CALL CLAMP(t, zero, one)
        CALL PUSHCONTROL3B(0)
      ELSE
        CALL DOT_PROD(c, d1, r)
        IF (e .LE. eps) THEN
          t = zero
          s = -(c/a)
          CALL PUSHREAL8(s)
          CALL CLAMP(s, zero, one)
          CALL PUSHCONTROL3B(1)
        ELSE
! General non-degenerate case
          CALL DOT_PROD(b, d1, d2)
          denom = a*e - b*b
          IF (denom .NE. zero) THEN
            s = (b*f-c*e)/denom
            CALL PUSHREAL8(s)
            CALL CLAMP(s, zero, one)
            CALL PUSHCONTROL1B(0)
          ELSE
            s = zero
            CALL PUSHCONTROL1B(1)
          END IF
          t = (b*s+f)/e
          IF (t .LT. zero) THEN
            t = zero
            CALL PUSHREAL8(s)
            s = -(c/a)
            CALL PUSHREAL8(s)
            CALL CLAMP(s, zero, one)
            CALL PUSHCONTROL3B(2)
          ELSE IF (t .GT. one) THEN
            t = one
            CALL PUSHREAL8(s)
            s = (b-c)/a
            CALL PUSHREAL8(s)
            CALL CLAMP(s, zero, one)
            CALL PUSHCONTROL3B(3)
          ELSE
            CALL PUSHCONTROL3B(4)
          END IF
        END IF
      END IF
      c1 = p1 + d1*s
      c2 = p2 + d2*t
      diff = c2 - c1
      dummy3 = diff
      diffb = 0.0_8
      dummy3b = 0.0_8
      CALL DOT_PROD_B(dsquared, dsquaredb, diff, diffb, dummy3, dummy3b)
      diffb = diffb + dummy3b
      c1b = 0.0_8
      c2b = 0.0_8
      c2b = diffb
      c1b = -diffb
      d2b = 0.0_8
      p2b = p2b + c2b
      d2b = t*c2b
      tb = SUM(d2*c2b)
      d1b = 0.0_8
      p1b = p1b + c1b
      d1b = s*c1b
      sb = SUM(d1*c1b)
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(t)
          CALL CLAMP_B(t, tb, zero, zerob, one, oneb)
          fb = tb/e
          eb = -(f*tb/e**2)
          zerob = zerob + sb
          rb = 0.0_8
          ab = 0.0_8
          GOTO 100
        ELSE
          CALL POPREAL8(s)
          CALL CLAMP_B(s, sb, zero, zerob, one, oneb)
          cb = -(sb/a)
          ab = c*sb/a**2
          zerob = zerob + tb
          eb = 0.0_8
          fb = 0.0_8
        END IF
      ELSE
        IF (branch .EQ. 2) THEN
          CALL POPREAL8(s)
          CALL CLAMP_B(s, sb, zero, zerob, one, oneb)
          CALL POPREAL8(s)
          cb = -(sb/a)
          ab = c*sb/a**2
          zerob = zerob + tb
          sb = 0.0_8
          tb = 0.0_8
          bb = 0.0_8
        ELSE IF (branch .EQ. 3) THEN
          CALL POPREAL8(s)
          CALL CLAMP_B(s, sb, zero, zerob, one, oneb)
          CALL POPREAL8(s)
          tempb1 = sb/a
          bb = tempb1
          cb = -tempb1
          ab = -((b-c)*tempb1/a)
          oneb = oneb + tb
          sb = 0.0_8
          tb = 0.0_8
        ELSE
          ab = 0.0_8
          bb = 0.0_8
          cb = 0.0_8
        END IF
        tempb0 = tb/e
        bb = bb + s*tempb0
        sb = sb + b*tempb0
        fb = tempb0
        eb = -((b*s+f)*tempb0/e)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(s)
          CALL CLAMP_B(s, sb, zero, zerob, one, oneb)
          tempb = sb/denom
          bb = bb + f*tempb
          fb = fb + b*tempb
          cb = cb - e*tempb
          eb = eb - c*tempb
          denomb = -((b*f-c*e)*tempb/denom)
        ELSE
          zerob = zerob + sb
          denomb = 0.0_8
        END IF
        ab = ab + e*denomb
        eb = eb + a*denomb
        bb = bb - 2*b*denomb
        CALL DOT_PROD_B(b, bb, d1, d1b, d2, d2b)
      END IF
      rb = 0.0_8
      CALL DOT_PROD_B(c, cb, d1, d1b, r, rb)
    END IF
 100 CALL DOT_PROD_B(f, fb, d2, d2b, r, rb)
    dummy3 = d2
    dummy3b = 0.0_8
    CALL DOT_PROD_B(e, eb, d2, d2b, dummy3, dummy3b)
    d2b = d2b + dummy3b
    dummy3 = d1
    dummy3b = 0.0_8
    CALL DOT_PROD_B(a, ab, d1, d1b, dummy3, dummy3b)
    d1b = d1b + dummy3b
    p1b = p1b + rb - d1b
    p2b = p2b - d2b - rb
    q2b = q2b + d2b
    q1b = q1b + d1b
    dsquaredb = 0.0_8
  END SUBROUTINE LINE_LINE_B
!  Differentiation of dot_prod in forward (tangent) mode:
!   variations   of useful results: d
!   with respect to varying inputs: v w
  SUBROUTINE DOT_PROD_D(d, dd, v, vd, w, wd)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: v, w
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: vd, wd
    REAL(kind=8), INTENT(OUT) :: d
    REAL(kind=8), INTENT(OUT) :: dd
    dd = vd(1)*w(1) + v(1)*wd(1) + vd(2)*w(2) + v(2)*wd(2) + vd(3)*w(3) &
&     + v(3)*wd(3)
    d = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
  END SUBROUTINE DOT_PROD_D
!  Differentiation of dot_prod in reverse (adjoint) mode:
!   gradient     of useful results: d v w
!   with respect to varying inputs: v w
  SUBROUTINE DOT_PROD_B(d, db, v, vb, w, wb)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: v, w
    REAL(kind=8), DIMENSION(3) :: vb, wb
    REAL(kind=8) :: d
    REAL(kind=8) :: db
    vb(1) = vb(1) + w(1)*db
    wb(1) = wb(1) + v(1)*db
    vb(2) = vb(2) + w(2)*db
    wb(2) = wb(2) + v(2)*db
    vb(3) = vb(3) + w(3)*db
    wb(3) = wb(3) + v(3)*db
  END SUBROUTINE DOT_PROD_B
!  Differentiation of clamp in forward (tangent) mode:
!   variations   of useful results: n
!   with respect to varying inputs: n min max
  SUBROUTINE CLAMP_D(n, nd, min, mind, max, maxd)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: min, max
    REAL(kind=8), INTENT(IN) :: mind, maxd
    REAL(kind=8), INTENT(INOUT) :: n
    REAL(kind=8), INTENT(INOUT) :: nd
    IF (n .LT. min) THEN
      nd = mind
      n = min
    END IF
    IF (n .GT. max) THEN
      nd = maxd
      n = max
    END IF
  END SUBROUTINE CLAMP_D
!  Differentiation of clamp in reverse (adjoint) mode:
!   gradient     of useful results: n min max
!   with respect to varying inputs: n min max
  SUBROUTINE CLAMP_B(n, nb, min, minb, max, maxb)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: min, max
    REAL(kind=8) :: minb, maxb
    REAL(kind=8), INTENT(INOUT) :: n
    REAL(kind=8) :: nb
    INTEGER :: branch
    IF (n .LT. min) THEN
      n = min
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (n .GT. max) THEN
      maxb = maxb + nb
      nb = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      minb = minb + nb
      nb = 0.0_8
    END IF
  END SUBROUTINE CLAMP_B
  SUBROUTINE LINE_LINE(p1, q1, p2, q2, dsquared)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: p1, q1, p2, q2
    REAL(kind=8), INTENT(OUT) :: dsquared
    REAL(kind=8), DIMENSION(3) :: d1, d2, r, diff, c1, c2, dummy3
    REAL(kind=8), PARAMETER :: eps=1e-12
    REAL(kind=8) :: a, b, c, e, f, s, t, denom
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    dummy3 = d1
    CALL DOT_PROD(a, d1, dummy3)
    dummy3 = d2
    CALL DOT_PROD(e, d2, dummy3)
    CALL DOT_PROD(f, d2, r)
    IF (a .LE. eps .AND. e .LE. eps) THEN
! both segments degenrate into points
      diff = q1 - p1
      dummy3 = diff
      CALL DOT_PROD(dsquared, diff, dummy3)
      RETURN
    ELSE
      IF (a .LE. eps) THEN
        s = zero
        t = f/e
        CALL CLAMP(t, zero, one)
      ELSE
        CALL DOT_PROD(c, d1, r)
        IF (e .LE. eps) THEN
          t = zero
          s = -(c/a)
          CALL CLAMP(s, zero, one)
        ELSE
! General non-degenerate case
          CALL DOT_PROD(b, d1, d2)
          denom = a*e - b*b
          IF (denom .NE. zero) THEN
            s = (b*f-c*e)/denom
            CALL CLAMP(s, zero, one)
          ELSE
            s = zero
          END IF
          t = (b*s+f)/e
          IF (t .LT. zero) THEN
            t = zero
            s = -(c/a)
            CALL CLAMP(s, zero, one)
          ELSE IF (t .GT. one) THEN
            t = one
            s = (b-c)/a
            CALL CLAMP(s, zero, one)
          END IF
        END IF
      END IF
      c1 = p1 + d1*s
      c2 = p2 + d2*t
      diff = c2 - c1
      dummy3 = diff
      CALL DOT_PROD(dsquared, diff, dummy3)
      RETURN
    END IF
  END SUBROUTINE LINE_LINE
  SUBROUTINE DOT_PROD(d, v, w)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: v, w
    REAL(kind=8), INTENT(OUT) :: d
    d = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
  END SUBROUTINE DOT_PROD
  SUBROUTINE CLAMP(n, min, max)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: min, max
    REAL(kind=8), INTENT(INOUT) :: n
    IF (n .LT. min) n = min
    IF (n .GT. max) n = max
  END SUBROUTINE CLAMP
END MODULE TRIANGLES_DB
