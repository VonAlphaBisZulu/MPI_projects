sup = [];
for p=prody_short
    sup = [sup ; p.sup];
end

slopermin = (([prody_short.YminPS_ramin]-[prody_short.YminPS_ra0])./[prody_short.ratpmin])';

slopermax = (([prody_short.YminPS_ramax]-[prody_short.YminPS_ra0])./[prody_short.ratpmax])';

a=(slopermax-slopermin);
z=zeros(length(prody_short),1);

YminPS = [prody_short.YminPA].*[prody_short.ratpPAmin]./10;

ratpmin = [z,[prody_short.ratpmin]',[prody_short.YminPS_ra0]',[prody_short.YminPS_ramin]'];
ratp    = [[prody_short.ratpmin]',[prody_short.ratpPAmin]',[prody_short.YminPS_ramin]',YminPS'];
ratpmax = [[prody_short.ratpPAmin]',[prody_short.ratpmax]',YminPS',[prody_short.YminPS_ramax]'];

clear metabolites
singleMCS = 1:length(prody_short);

names = cellstr(num2str((1:length(prody_short))'))';

plot([ratpmin(singleMCS,[1 2])'; ratp(singleMCS,2)'; ratpmax(singleMCS,2)'],[ratpmin(singleMCS,[3 4])' ;ratp(singleMCS,4)'; ratpmax(singleMCS,4)']);
hold on
plot([3.15 3.15],[0 2],'black');
hold off
legend(names);