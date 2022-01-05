function [Mi,nterms]=multivariateExponents(n,N)
%
%                          MATRIX EXPONENT
%
% This function build the indeces matrix for a multivariate polynomial of n
% variables and N degree (the two inputs). The vector nterms
% contain the number of possible indeces unique combinations (comibinations
% with repetitions) that are:
%
%      nterms=factorial(n+N-1)/factorial(N)/factorial(n-1)
%
% The matrix Mi contains all the unique combinations of exponents of the
% polynomial variables up at the desired order N.
%
% The used method is the following. At the order N, the following equation
% must be valid:
%
%                 e1 + e2 + ... + en = N   (1)
%
% where ej are the exponents for all the n variables that can be integers
% or zero. Thus, first the ordered unique solutions to this equations are
% found by menas of the function expomatr. This is a recursive function
% that generate a matrix em containing the solutions of the following:
%
%                 e1 + e2 + ... + eN = N   (2)
%                 e1>=e2>=e3>=....>=eN
%
% Then this matrix is either kept (N=n) or completed with zeros (N<n) or
% properly cutted in the desired way for (N>n) generating the matrix sols.
% For instance if N=3:
%
%                 em=[3 0 0; 2 1 0; 1 1 1]
%
% n=4 -->         sols=[3 0 0 0; 2 1 0 0; 1 1 1 0]
% n=3 -->         sols=em                                (3)
% n=2 -->         sols=[3 0; 2 1]
%
% This matrix sols represents the unique solutions of problem (1) with the
% coefficients satisfying the constaint: e1>=e2>=e3>=....>=en . From here,
% the constraint is taken away and so the solutions to equation (1) are the
% non-repeating permutations of each line of the matrix sols in (3) that
% are all collected in the matrix Mi. E.g. N=3, n=2
%
%                   Mi=[3 0; 0 3; 2 1; 1 2]
%
% The real bottleneck...
% of the algorithm is the permutation procedure. Still, a fast algorithm
% considering the repetitions in a smart way is not available. We can avoid
% the use of perms with several variables (it is inefficient with more than
% 10 variables) by using the combinations of indexes as done for the case
% of N<n.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix Initialization
Mi=[];
ii=N;
em=expomatr(ii);
% Split in the mentioned three cases
if ii==n
    sols=em;
else
    if ii<n
        [rr,cc]=size(em);
        sols=[em zeros(rr,n-cc)]; % Adjoining zeros
    else
        [rr,~]=size(em);
        % Cutting of the matrix: I need only the rows with no more than
        % n terms here
        vexct=(sum((em>0),2)<n+1);
        rindext=sort([1:rr]'.*vexct);
        rindex=rindext(rr-sum(vexct)+1:rr);
        cindex=1:n;
        sols=em(rindex,cindex);
    end
end
% This is the final solution of problem (2)
exppox=sortrows(sols,'descend');
[rr,~]=size(exppox);
% Now, each solution has to be permuted: in case ii<n there are a lot
% of zeros in the solutions. So, instead of using perms it is known
% that the permutation of the solution can be written as the
% combination of indexes of the variables taken nval at time where nval
% is the number of non zero terms in the solutions. Of course, this as
% to be repeated for every permutation of the non-zero nval terms of
% the solution. This procedure is not extended to the case ii>=n where
% the repeatitions can present with 1, 2 or other terms.
%
% Cycle over the solutions
for jj=1:rr
    if ii<n
        % Evaluating the non-zero terms of the current solution...
        nval=sum(exppox(jj,:)>0);
        % ...with their possible index combinations.
        combinind=nchoosek([1:n],nval);
        [rco,~]=size(combinind);
        % The selected solution is then permuted in the non-zero terms
        perval=sortrows(permsworep(exppox(jj,1:nval)),'descend');
        [rpe,~]=size(perval);
        % The Mi matrix is assembles by putting the indexes to each
        % term of the solution for every possible permutation.
        for kk=1:rpe
            Mit=zeros(rco,n);
            for ll=1:nval
                for mm=1:rco
                    Mit(mm,combinind(mm,ll))=perval(kk,ll);
                end
            end
            Mi=[Mi; Mit];
        end
    else
        Mi=[Mi; permsworep(exppox(jj,:))];
    end
end
[nterms,~]=size(Mi);
%     rth=factorial(n+ii-1)/factorial(ii)/factorial(n-1);
%     ntermsth=[ntermsth; rth];
Mi=sortrows(Mi,'descend');
end

function PnR=permsworep(v)

% PERMUTATIONS WITHOUT REPETITIONS
% Given a set of numbers in v, this code returns all the possible permuta-
% tions without repeating ones. Given that the length of vector v is n, all
% the possible permutations are n! in number and they are contained in the
% matrix P. This is computed with the matlab functions perms (that handles
% well up to n=10). If in v there are nr repeating terms i=1,...r each of
% them repeating ki times, the number of all possible permutations without
% repetitions is n!/k1!/k2!.../kr! and they are all saved in the output
% matrix PnR. These are obtained by properly sorting the matrix P coming
% out perms, evaluating the denominator k1!*k2!*...*kr! by examining the
% vector v and in the end taking 1 row of P every k1!*k2!*...*kr! ones.

n=length(v); % Dimension of the set v
P=perms(v); % P, matrix of all possible permutations

% Analysis of vector v, sorting and establishing repeating elements
vo=sort(v,'descend'); % Sorting
dvo=diff(vo);
% 2 Cases: if one element is zero, we have a repeating element, if not...
if sum(dvo==0)==0
    %PnR=sortrows(P,'descend'); % If there are not rep elemetents, PnR=P
    PnR=sortrows(P);
else
    % Here we need to count the number of reps per each rep element
    crep=1; % Rep. counter
    repv=[]; % Rep. vector, each line i is a rep term with repv(i) reps
    for ii=1:n-1
        if dvo(ii)==0
            crep=crep+1;
            % Particular case if we have the last element vo that is rep.
            if ii==n-1
                repv=[repv; crep];
            end
        else
            if crep>1
                repv=[repv; crep];
                crep=1;
            end
        end
    end
    step=prod(factorial(repv)); % Evaluation of k1!*k2!*...*kr!
    Po=sortrows(P,'descend'); % Sorting the rows of P
    PnR=Po(1:step:end,:); % Selecting every k1!*k2!*...*kr! row of Po
end
end

function E=expomatr(n)
sols=[];
if n==0
    E=[];
end
if n==1
    E=1;
end
if n==2
    E=[2 0; 1 1];
end
if n>2
    sols=[sols; n zeros(1,n-1)];
    for jj=n-1:-1:2
        dif=n-jj;
        Ejj=expomatr(dif);
        [rr,~]=size(Ejj);
        Ejjn=[];
        for kk=1:rr
            if Ejj(kk,1)<=jj
                Ejjn=[Ejjn; Ejj(kk,:)];
            end
        end
        Ejj=Ejjn;
        [rr,cc]=size(Ejj);
        sols=[sols; jj*ones(rr,1) Ejj zeros(rr,n-1-cc)];
    end
    sols=[sols; ones(1,n)];
    E=sortrows(sols);
end
end
