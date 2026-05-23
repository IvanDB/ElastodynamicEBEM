function[d, e, z, ierr]= exAlgolProcedure(n, d, e)
%     this subroutine is a translation of an algol procedure,
%     num. math. 12, 377-383(1968) by martin and wilkinson,
%     as modified in num. math. 15, 450(1970) by dubrulle.
%     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
%     this is a modified version of the 'eispack' routine imtql2.
%
%     this subroutine finds the eigenvalues and first components of the
%     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
%     method.
%
%     on input:
%
%        n is the order of the matrix;
%
%        d contains the diagonal elements of the input matrix;
%
%        e contains the subdiagonal elements of the input matrix
%          in its first n-1 positions.  e(n) is arbitrary;
%
%        z contains the first row of the identity matrix.
%
%      on output:
%
%        d contains the eigenvalues in ascending order.  if an
%          error exit is made, the eigenvalues are correct but
%          unordered for indices 1, 2, ..., ierr-1;
%
%        e has been destroyed;
%
%        z contains the first components of the orthonormal eigenvectors
%          of the symmetric tridiagonal matrix.  if an error exit is
%          made, z contains the eigenvectors associated with the stored
%          eigenvalues;
%
%        ierr is set to
%          zero       for normal return,
%          j          if the j-th eigenvalue has not been
%                     determined after 30 iterations.
%
%     ------------------------------------------------------------------
%
z=zeros(1,n); z(1)=1;

%
%
ierr = 0;
if  n ~= 1 
    %
    e(n) = 0;
    for l = 1 : n
        j = 0;
        %     :::::::::: look for small sub-diagonal element ::::::::::
        ivar =1;
        while ivar ==1
            for  m = l : n
                if m == n
                    break
                end
                if abs(e(m)) <= eps * (abs(d(m)) + abs(d(m+1)))
                    break
                end
                %          end  % if  m ~= n
            end % for m
            %
            p = d(l);
            if  m ~= l
                
                if j == 30
                    ierr=l;
                    return,
                end
                j = j + 1;
                %     :::::::::: form shift ::::::::::
                g = (d(l+1) - p) / (2 * e(l));
                r = sqrt(g*g+1);
                if g  < 0
                    ggg=-abs(r);
                else
                    ggg=abs(r);
                end
                g = d(m) - p + e(l) / (g + ggg);
                s = 1;
                c = 1;
                p = 0;
                mml = m - l;
                
                %     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
                
                
                for ii = 1: mml %      200 ii = 1, mml
                    i = m - ii;
                    f = s * e(i);
                    b = c * e(i);
                    if abs(f) >= abs(g)
                        c = g / f;
                        r = sqrt(c*c+1.0);
                        e(i+1) = f * r;
                        s = 1.0 / r;
                        c = c * s;
                    else
                        s = f / g;
                        r = sqrt(s*s+1.0);
                        e(i+1) = g * r;
                        c = 1.0 / r;
                        s = s * c;
                    end % if
                    g = d(i+1) - p;
                    r = (d(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d(i+1) = g + p;
                    g = c * r - b;
                    %     :::::::::: form first component of vector ::::::::::
                    f = z(i+1);
                    z(i+1) = s * z(i) + c * f;
                    z(i) = c * z(i) - s * f;
                end % for
                %
                d(l) = d(l) - p;
                e(l) = g;
                e(m) = 0.0;
            else
                ivar =2;
                
            end % if ivar
            %   end % forse questo non ci vuuole
        end % while
    end % for        240 continue
    %
    %     :::::::::: order eigenvalues and eigenvectors ::::::::::
    for  ii = 2 : n
        i = ii - 1;
        k = i;
        p = d(i);
        %
        for  j = ii : n
            if d(j) <= p
                k = j;
                p = d(j);
            end % if
        end % for   260    continue
        %
        if k ~= i
            d(k) = d(i);
            d(i) = p;
            p = z(i);
            z(i) = z(k);
            z(k) = p;
        end % if
    end % for
    %
    %     :::::::::: set error -- no convergence to an
    %                eigenvalue after 30 iterations ::::::::::
end % if      1001 return
%c     :::::::::: last card of gausq2 ::::::::::
