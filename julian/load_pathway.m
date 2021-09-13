function [prod_id,species,reactions] = load_pathway(product)
    switch product
        case  1 % ethanol
			prod_id = 'EX_etoh_e';
			species = [];
			reactions = [];
        case  2 % lysine
			prod_id      = 'EX2_lys__L_e';
			species      = [];
			reactions(1) = struct('reac_id','EX2_lys__L_e','equation','1 lys__L_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  3 % glutamate
			prod_id      = 'EX2_glu__L_e';
			species      = [];
			reactions(1) = struct('reac_id','EX2_glu__L_e','equation','1 glu__L_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  4 % isobutanol
			prod_id      = 'EX_ibutoh_e';
			species(1)   = struct('spec_id','2mppal_c','spec_name','2-Methylpropanal','fbc_chemicalFormula','C4H8O','fbc_charge',0);
			species(2)   = struct('spec_id','ibutoh_c','spec_name','Isobutyl alcohol','fbc_chemicalFormula','C4H10O','fbc_charge',0);
			species(3)   = struct('spec_id','ibutoh_p','spec_name','Isobutyl alcohol','fbc_chemicalFormula','C4H10O','fbc_charge',0);
			species(4)   = struct('spec_id','ibutoh_e','spec_name','Isobutyl alcohol','fbc_chemicalFormula','C4H10O','fbc_charge',0);
			reactions(1) = struct('reac_id','3MOBDC','equation','1 3mob_c + 1 h_c = 1 co2_c + 1 2mppal_c','lb',0,'ub',1000,'fbc_geneProductAssociation','kivD');
			reactions(2) = struct('reac_id','ALCD23xi','equation','1 h_c + 1 nadh_c + 1 2mppal_c = 1 nad_c + 1 ibutoh_c','lb',0,'ub',1000,'fbc_geneProductAssociation','adh2');
			reactions(3) = struct('reac_id','IBUTOHtpp','equation','1 ibutoh_c = 1 ibutoh_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(4) = struct('reac_id','IBUTOHtex','equation','1 ibutoh_p = 1 ibutoh_e','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(5) = struct('reac_id','EX_ibutoh_e','equation','1 ibutoh_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  5 % 1,4-BDO
			prod_id      = 'EX_14bdo_e';
			species(1)   = struct('spec_id','4hb_c','spec_name','4-Hydroxybutanoate','fbc_chemicalFormula','C4H7O3','fbc_charge',-1);
			species(2)   = struct('spec_id','4hbcoa_c','spec_name','4-Hydroxybutyryl-CoA','fbc_chemicalFormula','C25H38N7O18P3S','fbc_charge',-4);
			species(3)   = struct('spec_id','4hbal_c','spec_name','4-Hydroxybutanal','fbc_chemicalFormula','C4H8O2','fbc_charge',0);
			species(4)   = struct('spec_id','14bdo_c','spec_name','Butane-1,4-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			species(5)   = struct('spec_id','14bdo_p','spec_name','Butane-1,4-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			species(6)   = struct('spec_id','14bdo_e','spec_name','Butane-1,4-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			reactions(1) = struct('reac_id','SSCOARx','equation','1 h_c + 1 nadph_c + 1 succoa_c = 1 coa_c + 1 nadp_c + 1 sucsal_c','lb',0,'ub',1000,'fbc_geneProductAssociation','sucD_Pg');
			reactions(2) = struct('reac_id','AKGDC','equation','1 akg_c + 1 h_c = 1 co2_c + 1 sucsal_c','lb',0,'ub',1000,'fbc_geneProductAssociation','sucA_Mb');
			reactions(3) = struct('reac_id','4HBD','equation','1 h_c + 1 nadh_c + 1 sucsal_c = 1 4hb_c + 1 nad_c','lb',0,'ub',1000,'fbc_geneProductAssociation','4hbD_Pg');
			reactions(4) = struct('reac_id','4HBCT','equation','1 4hb_c + 1 accoa_c = 1 4hbcoa_c + 1 ac_c','lb',0,'ub',1000,'fbc_geneProductAssociation','cat2');
			reactions(5) = struct('reac_id','4HBDH','equation','1 4hbcoa_c + 1 h_c + 1 nadh_c = 1 4hbal_c + 1 coa_c + 1 nad_c','lb',0,'ub',1000,'fbc_geneProductAssociation','adh2');
			reactions(6) = struct('reac_id','4HBDx','equation','1 4hbal_c + 1 h_c + 1 nadh_c = 1 14bdo_c + 1 nad_c','lb',0,'ub',1000,'fbc_geneProductAssociation','adhE2');
			reactions(7) = struct('reac_id','14BDOtpp','equation','1 14bdo_c = 1 14bdo_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(8) = struct('reac_id','14BDOtex','equation','1 14bdo_p = 1 14bdo_e','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(9) = struct('reac_id','EX_14bdo_e','equation','1 14bdo_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  6 % 2,3-BDO
			prod_id      = 'EX_23bdo_e';
			species(1)   = struct('spec_id','actn_c','spec_name','3-Hydroxybutan-2-one','fbc_chemicalFormula','C4H8O2','fbc_charge',0);
			species(2)   = struct('spec_id','23bdo_c','spec_name','Butane-2,3-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			species(3)   = struct('spec_id','23bdo_p','spec_name','Butane-2,3-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			species(4)   = struct('spec_id','23bdo_e','spec_name','Butane-2,3-diol','fbc_chemicalFormula','C4H10O2','fbc_charge',0);
			reactions(1) = struct('reac_id','ACLDC','equation','1 alac__S_c + 1 h_c = 1 co2_c + 1 actn_c','lb',0,'ub',1000,'fbc_geneProductAssociation','aldA and aldC');
			reactions(2) = struct('reac_id','BTDD','equation','1 h_c + 1 nadh_c + 1 actn_c = 1 nad_c + 1 23bdo_c','lb',0,'ub',1000,'fbc_geneProductAssociation','butB');
			reactions(3) = struct('reac_id','23BDOtpp','equation','1 23bdo_c = 1 23bdo_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(4) = struct('reac_id','23BDOtex','equation','1 23bdo_p = 1 23bdo_e','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(5) = struct('reac_id','EX_23bdo_e','equation','1 23bdo_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  7 % itaconic acid
			prod_id = 'EX_itacon_e';
			species(1)   = struct('spec_id','itacon_c','spec_name','2-Methylenesuccinate','fbc_chemicalFormula','C5H4O4','fbc_charge',-2);
			species(2)   = struct('spec_id','itacon_p','spec_name','2-Methylenesuccinate','fbc_chemicalFormula','C5H4O4','fbc_charge',-2);
			species(3)   = struct('spec_id','itacon_e','spec_name','2-Methylenesuccinate','fbc_chemicalFormula','C5H4O4','fbc_charge',-2);
			reactions(1) = struct('reac_id','CADC','equation','1 h_c + 1 acon_C_c = 1 itacon_c + 1 co2_c','lb',0,'ub',1000,'fbc_geneProductAssociation','cadA');
			reactions(2) = struct('reac_id','ITACONtpp','equation','1 itacon_c = 1 itacon_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(3) = struct('reac_id','ITACONtex','equation','1 itacon_p = 1 itacon_e','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(4) = struct('reac_id','EX_itacon_e','equation','1 itacon_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case  8 % isoprene
			prod_id      = 'EX_is_e';
			species(1)   = struct('spec_id','hmgcoa_c','spec_name','3-Hydroxy-3-methyl-glutaryl-CoA','fbc_chemicalFormula','C27H39N7O20P3S','fbc_charge',-5);
			species(2)   = struct('spec_id','mev_c','spec_name','(3R)-3,5-Dihydroxy-3-methylpentanoate','fbc_chemicalFormula','C6H11O4','fbc_charge',-1);
			species(3)   = struct('spec_id','5pmev_c','spec_name','(3R)-3-Hydroxy-3-methyl-5-(phosphonooxy)pentanoate','fbc_chemicalFormula','C6H10O7P','fbc_charge',-3);
			species(4)   = struct('spec_id','5dpmev_c','spec_name','(3R)-3-Hydroxy-5-{[hydroxy(phosphonooxy)phosphoryl]oxy}-3-methylpentanoate','fbc_chemicalFormula','C6H10O10P2','fbc_charge',-4);
			species(5)   = struct('spec_id','is_c','spec_name','Isoprene','fbc_chemicalFormula','C5H8','fbc_charge',0);
			species(6)   = struct('spec_id','is_p','spec_name','Isoprene','fbc_chemicalFormula','C5H8','fbc_charge',0);
			species(7)   = struct('spec_id','is_e','spec_name','Isoprene','fbc_chemicalFormula','C5H8','fbc_charge',0);
			reactions(1) = struct('reac_id','HMGCOAS','equation','1 aacoa_c + 1 accoa_c + 1 h2o_c = 1 coa_c + 1 h_c + 1 hmgcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','mvaS');
			reactions(2) = struct('reac_id','HMGCOAR','equation','2 h_c + 2 nadph_c + 1 hmgcoa_c = 1 coa_c + 2 nadp_c + 1 mev_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','mvaE');
			reactions(3) = struct('reac_id','MEVK1','equation','1 atp_c + 1 mev_c = 1 adp_c + 1 h_c + 1 5pmev_c','lb',0,'ub',1000,'fbc_geneProductAssociation','mvk');
			reactions(4) = struct('reac_id','PMEVK','equation','1 atp_c + 1 5pmev_c = 1 adp_c + 1 5dpmev_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pmk');
			reactions(5) = struct('reac_id','DPMVD','equation','1 atp_c + 1 5dpmev_c = 1 adp_c + 1 co2_c + 1 ipdp_c + 1 pi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','mvd');
			reactions(6) = struct('reac_id','ISS','equation','1 dmpp_c = 1 is_c + 1 ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','IspS');
			reactions(7) = struct('reac_id','IStpp','equation','1 is_c = 1 is_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(8) = struct('reac_id','IStex','equation','1 is_p = 1 is_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(9) = struct('reac_id','EX_is_e','equation','1 is_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case 9 % butane
			prod_id      = 'EX_butane_e';
            species(1)   = struct('spec_id','ppoh_c','spec_name','','fbc_chemicalFormula','C3H8O1','fbc_charge',0);
			species(2)   = struct('spec_id','val_c','spec_name','','fbc_chemicalFormula','C5H9O2','fbc_charge',-1);
			species(3)   = struct('spec_id','ptcoa_c','spec_name','','fbc_chemicalFormula','C26H40N7O17P3S','fbc_charge',-4);
			species(4)   = struct('spec_id','cis23dhbtcoa_c','spec_name','','fbc_chemicalFormula','C25H36N7O17P3S','fbc_charge',-4);
			species(5)   = struct('spec_id','cis2mctcoa_c','spec_name','','fbc_chemicalFormula','C26H38N7O17P3S','fbc_charge',-4);
			species(6)   = struct('spec_id','3hbcoa__R_c','spec_name','','fbc_chemicalFormula','C25H38N7O18P3S','fbc_charge',-4);
			species(7)   = struct('spec_id','3hptcoa_c','spec_name','','fbc_chemicalFormula','C26H40N7O18P3S','fbc_charge',-4);
			species(8)   = struct('spec_id','3oxptcoa_c','spec_name','','fbc_chemicalFormula','C26H38N7O18P3S','fbc_charge',-4);
            species(9)   = struct('spec_id','ptal_c','spec_name','','fbc_chemicalFormula','C5H10O','fbc_charge',0);
            species(10)   = struct('spec_id','butane_c','spec_name','','fbc_chemicalFormula','C4H10','fbc_charge',0);
			species(11)   = struct('spec_id','butane_p','spec_name','','fbc_chemicalFormula','C4H10','fbc_charge',0);
			species(12)   = struct('spec_id','butane_e','spec_name','','fbc_chemicalFormula','C4H10','fbc_charge',0);
			reactions(1) = struct('reac_id','ALCD3ir','equation','1 h_c + 1 ppal_c + 1 nadh_c = 1 nad_c + 1 ppoh_c','lb',0,'ub',0,'fbc_geneProductAssociation','E2348C_RS11505');
			reactions(2) = struct('reac_id','MMSAD2','equation','1 nad_c + 1 ppal_c + 1 coa_c = 1 h_c + 1 ppcoa_c + 1 nadh_c','lb',0,'ub',0,'fbc_geneProductAssociation','E2348C_RS11500');
			reactions(3) = struct('reac_id','ADO50','equation','1 h2o_c + 1 ptal_c = 1 h_c + 1 for_c + 1 butane_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','ado');
			reactions(4) = struct('reac_id','CARNi5','equation','1 h_c + 1 atp_c + 1 nadph_c + 1 val_c = 1 ppi_c + 1 nadp_c + 1 amp_c + 1 ptal_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','car');
			reactions(5) = struct('reac_id','THSTRS5','equation','1 h2o_c + 1 ptcoa_c = 1 h_c + 1 coa_c + 1 val_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','tesB');
			reactions(6) = struct('reac_id','TERTD4','equation','1 h_c + 1 nadh_c + 1 cis23dhbtcoa_c = 1 nad_c + 1 btcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','ter');
			reactions(7) = struct('reac_id','TERTD5','equation','1 h_c + 1 nadh_c + 1 cis2mctcoa_c = 1 nad_c + 1 ptcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','ter');
			reactions(8) = struct('reac_id','PHAJ4B4','equation','1 3hbcoa__R_c = 1 h2o_c + 1 cis23dhbtcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','phaJ4b');
			reactions(9) = struct('reac_id','PHAJ4B5','equation','1 3hptcoa_c = 1 h2o_c + 1 cis2mctcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','phaJ4b');
			reactions(10) = struct('reac_id','PHAB5','equation','1 h_c + 1 nadph_c + 1 3oxptcoa_c = 1 nadp_c + 1 3hptcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','phaB');
			reactions(11) = struct('reac_id','PHAB4','equation','1 h_c + 1 aacoa_c + 1 nadph_c = 1 nadp_c + 1 3hbcoa__R_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','phaB');
			reactions(12) = struct('reac_id','BTKB5','equation','1 accoa_c + 1 ppcoa_c = 1 coa_c + 1 3oxptcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','bktB');
			reactions(13) = struct('reac_id','BUTANEtp','equation','1 butane_c = 1 butane_p','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(14) = struct('reac_id','BUTANEtex','equation','1 butane_p = 1 butane_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(15) = struct('reac_id','EX_butane_e','equation','1 butane_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case 10 % methacrylic acid
			prod_id      = 'EX_2mp2e_e';
			species(1)   = struct('spec_id','2mppal_c','spec_name','2-Methylpropanal','fbc_chemicalFormula','C4H8O','fbc_charge',0);
			species(2)   = struct('spec_id','2mpp_c','spec_name','2-Methylpropanoate','fbc_chemicalFormula','C4H7O2','fbc_charge',-1);
			species(3)   = struct('spec_id','2mpcoa_c','spec_name','2-Methylpropanoly-CoA','fbc_chemicalFormula','C25H38N7O17P3S','fbc_charge',-4);
			species(4)   = struct('spec_id','2mp2coa_c','spec_name','2-Methylprop-2-enoyl-CoA','fbc_chemicalFormula','C25H36N7O17P3S','fbc_charge',-4);
			species(5)   = struct('spec_id','2mp2e_c','spec_name','2-Methylprop-2-enoic acid','fbc_chemicalFormula','C4H5O2','fbc_charge',-1);
			species(6)   = struct('spec_id','2mp2e_p','spec_name','2-Methylprop-2-enoic acid','fbc_chemicalFormula','C4H5O2','fbc_charge',-1);
			species(7)   = struct('spec_id','2mp2e_e','spec_name','2-Methylprop-2-enoic acid','fbc_chemicalFormula','C4H5O2','fbc_charge',-1);
			reactions(1) = struct('reac_id','3MOBDC','equation','1 3mob_c + 1 h_c = 1 co2_c + 1 2mppal_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pdc_1_5_6');
			reactions(2) = struct('reac_id','2MPD','equation','1 2mppal_c + 1 nad_c + h2o_c = 1 2mpp_c + 2 h_c + 1 nadh_c','lb',0,'ub',1000,'fbc_geneProductAssociation','padA');
			reactions(3) = struct('reac_id','2MPCT','equation','1 2mpp_c + 1 accoa_c = 1 2mpcoa_c + 1 ac_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pct');
			reactions(4) = struct('reac_id','BCDH','equation','1 3mob_c + 1 nad_c + coa_c = 1 2mpcoa_c + 1 co2_c + 1 nadh_c','lb',0,'ub',1000,'fbc_geneProductAssociation','bcdh');
			reactions(5) = struct('reac_id','SCDH','equation','1 fad_c + 2mpcoa_c = fadh2_c + 2mp2coa_c','lb',0,'ub',1000,'fbc_geneProductAssociation','acdH');
			reactions(6) = struct('reac_id','MTE','equation','1 h2o_c + 1 2mp2coa_c = 1 2mp2e_c + 1 coa_c + 1 h_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fcbC');
			reactions(7) = struct('reac_id','2MP2Etpp','equation','1 2mp2e_c = 1 2mp2e_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(8) = struct('reac_id','2MP2Etex','equation','1 2mp2e_p = 1 2mp2e_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(9) = struct('reac_id','EX_2mp2e_e','equation','1 2mp2e_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case 11 % resveratrol
			prod_id      = 'EX_revtl_e';
			species(1)   = struct('spec_id','cou_c','spec_name','4-coumarate','fbc_chemicalFormula','C9H7O3','fbc_charge',-1);
			species(2)   = struct('spec_id','coucoa_c','spec_name','4-coumaroyl-CoA','fbc_chemicalFormula','C30H38N7O18P3S','fbc_charge',-4);
			species(3)   = struct('spec_id','revtl_c','spec_name','trans-resveratrol','fbc_chemicalFormula','C14H12O3','fbc_charge',0);
			species(4)   = struct('spec_id','revtl_p','spec_name','trans-resveratrol (periplasm)','fbc_chemicalFormula','C14H12O3','fbc_charge',0);
			species(5)   = struct('spec_id','revtl_e','spec_name','trans-resveratrol (extracellular)','fbc_chemicalFormula','C14H12O3','fbc_charge',0);
			reactions(1) = struct('reac_id','TAL','equation','tyr__L_c = nh4_c + cou_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pal');
			reactions(2) = struct('reac_id','CCL','equation','1 cou_c + 1 atp_c + 1 coa_c = 1 coucoa_c + 1 amp_c + 1 ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','couB');
			reactions(3) = struct('reac_id','STS','equation','3 malcoa_c + coucoa_c + 3 h_c = 4 coa_c + revtl_c + 4 co2_c','lb',0,'ub',1000,'fbc_geneProductAssociation','sts');
			reactions(4) = struct('reac_id','REVTLtpp','equation','1 revtl_c = 1 revtl_p','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(5) = struct('reac_id','REVTLtex','equation','1 revtl_p = 1 revtl_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(6) = struct('reac_id','EX_revtl_e','equation','1 revtl_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(7) = struct('reac_id','XPK','equation','1 xu5p__D_c + 1 pi_c = 1 actp_c + 1 g3p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
			reactions(8) = struct('reac_id','FPK','equation','1 f6p_c + 1 pi_c = 1 actp_c + 1 e4p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
        case 12 % bisabolene
			prod_id      = 'EX_bsb_e';
			species(1)   = struct('spec_id','hmgcoa_c','spec_name','3-Hydroxy-3-methyl-glutaryl-CoA','fbc_chemicalFormula','C27H39N7O20P3S','fbc_charge',-5);
			species(2)   = struct('spec_id','mev_c','spec_name','(3R)-3,5-Dihydroxy-3-methylpentanoate','fbc_chemicalFormula','C6H11O4','fbc_charge',-1);
			species(3)   = struct('spec_id','5pmev_c','spec_name','(3R)-3-Hydroxy-3-methyl-5-(phosphonooxy)pentanoate','fbc_chemicalFormula','C6H10O7P','fbc_charge',-3);
			species(4)   = struct('spec_id','5dpmev_c','spec_name','(3R)-3-Hydroxy-5-{[hydroxy(phosphonooxy)phosphoryl]oxy}-3-methylpentanoate','fbc_chemicalFormula','C6H10O10P2','fbc_charge',-4);
			species(5)   = struct('spec_id','bsb_c','spec_name','bisabolene','fbc_chemicalFormula','C15H24','fbc_charge',0);
			species(6)   = struct('spec_id','bsb_p','spec_name','bisabolene (periplasm)','fbc_chemicalFormula','C15H24','fbc_charge',0);
			species(7)   = struct('spec_id','bsb_e','spec_name','bisabolene (extracellular)','fbc_chemicalFormula','C15H24','fbc_charge',0);
			reactions(1) = struct('reac_id','HMGCOAS','equation','1 aacoa_c + 1 accoa_c + 1 h2o_c = 1 coa_c + 1 h_c + 1 hmgcoa_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','mep');
			reactions(2) = struct('reac_id','HMGCOAR','equation','2 h_c + 2 nadph_c + 1 hmgcoa_c = 1 coa_c + 2 nadp_c + 1 mev_c','lb',-1000,'ub',1000,'fbc_geneProductAssociation','mep');
			reactions(3) = struct('reac_id','MEVK1','equation','1 atp_c + 1 mev_c = 1 adp_c + 1 h_c + 1 5pmev_c','lb',0,'ub',1000,'fbc_geneProductAssociation','mep');
			reactions(4) = struct('reac_id','PMEVK','equation','1 atp_c + 1 5pmev_c = 1 adp_c + 1 5dpmev_c','lb',0,'ub',1000,'fbc_geneProductAssociation','mep');
			reactions(5) = struct('reac_id','DPMVD','equation','1 atp_c + 1 5dpmev_c = 1 adp_c + 1 co2_c + 1 ipdp_c + 1 pi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','mep');
			reactions(6) = struct('reac_id','FPPS','equation','2 ipdp_c + 1 dmpp_c = frdp_c + 2 ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','ispA');
			reactions(7) = struct('reac_id','BSBS','equation','1 frdp_c = 1 bsb_c + 1 ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','ag1');
			reactions(8) = struct('reac_id','BSBtpp','equation','1 bsb_c = 1 bsb_p','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(9) = struct('reac_id','BSBtex','equation','1 bsb_p = 1 bsb_e','lb',-1000,'ub',1000,'fbc_geneProductAssociation','');
			reactions(10) = struct('reac_id','EX_bsb_e','equation','1 bsb_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
% 			reactions(11) = struct('reac_id','XPK','equation','1 xu5p__D_c + 1 pi_c = 1 actp_c + 1 g3p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
% 			reactions(12) = struct('reac_id','FPK','equation','1 f6p_c + 1 pi_c = 1 actp_c + 1 e4p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
        case  13 % tryptophane
			prod_id      = 'EX2_trp__L_e';
			species      = [];
			reactions(1) = struct('reac_id','EX2_trp__L_e','equation','1 trp__L_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
        case 14 % octyl acetate
			prod_id      = 'EX_ocac_e';
			species(1)   = struct('spec_id','octal_c','spec_name','octanal','fbc_chemicalFormula','C8H16O','fbc_charge',0);
			species(2)   = struct('spec_id','octoh_c','spec_name','octanol','fbc_chemicalFormula','C8H18O','fbc_charge',0);
			species(3)   = struct('spec_id','ocac_c','spec_name','octyl acetate','fbc_chemicalFormula','C10H20O2','fbc_charge',0);
			species(4)   = struct('spec_id','ocac_p','spec_name','octyl acetate (periplasm)','fbc_chemicalFormula','C10H20O2','fbc_charge',0);
			species(5)   = struct('spec_id','ocac_e','spec_name','octyl acetate (extracellular)','fbc_chemicalFormula','C10H20O2','fbc_charge',0);
			reactions(1) = struct('reac_id','CAR','equation','1 octa_c + nadph_c + h_c + atp_c = octal_c + nadp_c + amp_c + ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','carholo');
			reactions(2) = struct('reac_id','AHR','equation','1 octal_c + nadph_c + h_c = 1 octoh_c + nadp_c','lb',0,'ub',1000,'fbc_geneProductAssociation','oahr');
			reactions(3) = struct('reac_id','AAT','equation','1 octoh_c + 1 accoa_c = 1 ocac_c + 1 coa_c','lb',0,'ub',1000,'fbc_geneProductAssociation','scatf1');
			reactions(4) = struct('reac_id','OCACtpp','equation','1 ocac_c = 1 ocac_p','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(5) = struct('reac_id','OCACtex','equation','1 ocac_p = 1 ocac_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(6) = struct('reac_id','EX_ocac_e','equation','1 ocac_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
% 			reactions(7) = struct('reac_id','XPK','equation','1 xu5p__D_c + 1 pi_c = 1 actp_c + 1 g3p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
% 			reactions(8) = struct('reac_id','FPK','equation','1 f6p_c + 1 pi_c = 1 actp_c + 1 e4p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
        case 15 % 4-hydroxycoumarin
			prod_id      = 'EX_4hc_e';
			species(1)   = struct('spec_id','salc_c','spec_name','salicylate','fbc_chemicalFormula','C7H5O3','fbc_charge',-1);
			species(2)   = struct('spec_id','salccoa_c','spec_name','salicyl-CoA, 2-Hydroxybenzoyl-CoA','fbc_chemicalFormula','C28H36N7O18P3S','fbc_charge',-4);
			species(3)   = struct('spec_id','4hc_c','spec_name','4-hydroxycoumarin','fbc_chemicalFormula','C9H6O3','fbc_charge',0);
			species(4)   = struct('spec_id','4hc_p','spec_name','4-hydroxycoumarin (periplasm)','fbc_chemicalFormula','C9H6O3','fbc_charge',0);
			species(5)   = struct('spec_id','4hc_e','spec_name','4-hydroxycoumarin (extracellular)','fbc_chemicalFormula','C9H6O3','fbc_charge',0);
			reactions(1) = struct('reac_id','IPL','equation','ichor_c = pyr_c + salc_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pchB');
			reactions(2) = struct('reac_id','SCL','equation','salc_c + coa_c + atp_c = salccoa_c + amp_c + ppi_c','lb',0,'ub',1000,'fbc_geneProductAssociation','sdgA');
			reactions(3) = struct('reac_id','BIS','equation','salccoa_c + h_c + malcoa_c = 4hc_c + co2_c + 2 coa_c','lb',0,'ub',1000,'fbc_geneProductAssociation','pqsD');
			reactions(4) = struct('reac_id','4HCtpp','equation','1 4hc_c = 1 4hc_p','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(5) = struct('reac_id','4HCtex','equation','1 4hc_p = 1 4hc_e','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(6) = struct('reac_id','EX_4hc_e','equation','1 4hc_e =','lb',0,'ub',1000,'fbc_geneProductAssociation','');
			reactions(7) = struct('reac_id','XPK','equation','1 xu5p__D_c + 1 pi_c = 1 actp_c + 1 g3p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
			reactions(8) = struct('reac_id','FPK','equation','1 f6p_c + 1 pi_c = 1 actp_c + 1 e4p_c + 1 h2o_c','lb',0,'ub',1000,'fbc_geneProductAssociation','fxpk');
    end
end

