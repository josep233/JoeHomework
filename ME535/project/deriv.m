% Find expressions for lift and moment due to Theodorsen's solution.
%
syms pi rho b w h U al a C k psiB psiT x

% Lift force from Dowell et al.
Leq=-pi*rho*b^2*(-w^2*h+1i*w*U*al+b*a*w^2*al)...
    -2*pi*rho*U*b*C*(1i*w*h+U*al+b*(0.5-a)*1i*w*al);

Leq_h = simplify(subs(Leq,al,0));
Leq_al = simplify(subs(Leq,h,0));

% Leq_adj = [-psiB, -x*psiT] * [Leq_h; Leq_al];
Leq_adj = -psiB * Leq_h - x * psiT * Leq_al;

Leq = Leq_adj;

% pretty(Leq)
Leq=subs(Leq,U,b*w/k);
Leq=simplify(Leq/w^2); % Clear out all w^2 terms, since they appear in every expression.
Leq=simplify(Leq/(pi*rho*b^2)); % clear out pi*rho*b^2 from all terms
Lh=simplify(subs(Leq,al,0));
Lal=simplify(subs(Leq,h,0));

% Aerodynamic Moment from Dowell et al.
Meq=pi*rho*b^2*(-b*a*w^2*h-1i*w*U*b*(0.5-a)*al+b^2*(1/8-a^2)*w^2*al)+...
    2*pi*rho*U*b^2*(0.5-a)*C*(1i*w*h+U*al+b*(0.5-a)*1i*w*al);

Meq_h = simplify(subs(Meq,al,0));
Meq_al = simplify(subs(Meq,h,0));

% Meq_adj = [0, psiT] * [Meq_h; Meq_al];
Meq_adj = Meq_al * psiT;

Meq = Meq_adj;

Meq=subs(Meq,U,b*w/k);
Meq=simplify(Meq/w^2); % all terms need omega^2
Meq=simplify(Meq/(pi*rho*b^2)); % all terms need pi*rho*b^2
Mh=simplify(subs(Meq,al,0));
Mal=simplify(subs(Meq,h,0))
    simplify(Mal,'All',true)

Lh