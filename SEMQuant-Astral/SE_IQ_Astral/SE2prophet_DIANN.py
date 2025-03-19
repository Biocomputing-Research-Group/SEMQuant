import sys
class PSM:
    def __init__(self, filename, file, scan, ParentCharge, MeasuredParentMass,CalculatedParentMass,Massdiff, score,  IdentifiedPeptide,OriginalPeptide,
                 Proteins,ProteinCount):
        self.filename=filename
        self.file = file
        self.scan = scan
        self.ParentCharge = ParentCharge
        self.rank = 0
        self.MeasuredParentMass=MeasuredParentMass
        self.CalculatedParentMass=CalculatedParentMass
        self.Massdiff = Massdiff
        self.score = score
        self.IdentifiedPeptide = IdentifiedPeptide
        self.OriginalPeptide=OriginalPeptide
        self.Proteins = Proteins
        self.ProteinCount=ProteinCount

class Spectrum:
    def __init__(self,spectrumID,scan,precursor_neutral_mass,charge,PSM):
        self.spectrumID=spectrumID
        self.scan=scan
        self.precursor_neutral_mass=precursor_neutral_mass
        self.charge=charge
        self.PSMs=[PSM]

def sipros_parse(filename):
    spectrum_d=dict()
    with open(filename) as f:
        for line_id,line in enumerate(f):
            if line_id<60:
                continue
            s=line.strip().split('\t')
            filename=s[0]
            file=s[0].replace('.ms2','')
            scan=s[1]
            charge=s[2]
            MeasuredParentMass=float(s[3])
            CalculatedParentMass=float(s[4])
            massdiff=str(abs(MeasuredParentMass-CalculatedParentMass))
            score=s[10]
            IdentifiedPeptide=s[13].replace('[','')
            IdentifiedPeptide=IdentifiedPeptide.replace(']','')
            OriginalPeptide=s[14].replace('[','')
            OriginalPeptide=OriginalPeptide.replace(']','')
            Proteins=s[15].replace('{','')
            Proteins=Proteins.replace('}','')
            Proteins=Proteins.split(',')
            idx=file+'.'+scan+'.'+scan+'.'+charge
            PSM_ins = PSM(filename, file, scan, charge, MeasuredParentMass, CalculatedParentMass, massdiff, score,
                          IdentifiedPeptide, OriginalPeptide, Proteins, len(Proteins))
            if idx in spectrum_d:
                PSM_ins.rank=len(spectrum_d[idx].PSMs)+1
                spectrum_d[idx].PSMs.append(PSM)
            else:
                PSM_ins.rank=1
                Spectrum_ins=Spectrum(idx,scan,MeasuredParentMass,charge,PSM_ins)
                spectrum_d[idx]=Spectrum_ins



    return spectrum_d

def read_tab(tab_file):
    scan2RT={}
    with open(tab_file) as f:
        for line_id, line in enumerate(f):
            if line_id>0:
                s=line.strip().split('\t')
                scan2RT[s[1]]=str(float(s[-2])*60) # transfer from mins to secs
    return scan2RT



if __name__ == "__main__":
    spectrum_d=sipros_parse(sys.argv[1])
    scan2RT=read_tab(sys.argv[2]+'/'+sys.argv[3]+'.SE.tab')
    with open('/home/UNT/fs0199/SEMQuant/SE2IQ/ProteinProphetHeader_DIANN.txt') as f:
        header=f.read()

    with open(sys.argv[2]+'/'+sys.argv[3]+'.pep.xml','w') as f:
        f.write(header)
        for spectrum in spectrum_d:
            spectrumItem=spectrum_d[spectrum]
            spectrumItemString='<spectrum_query spectrum="{0}" start_scan="{1}" end_scan="{1}" precursor_neutral_mass="{2}" assumed_charge="{3}" index="0" retention_time_sec="{4}">\n'.format(spectrumItem.spectrumID,spectrumItem.scan,spectrumItem.precursor_neutral_mass,spectrumItem.charge,scan2RT[spectrumItem.scan])+'<search_result>\n'
            f.write(spectrumItemString)
            PSM_element=spectrumItem.PSMs
            for p in PSM_element:
                tests=p.OriginalPeptide[0:-1]
                num_missing_cleavages=tests.count('K')+tests.count('R')
                psmItemString='<search_hit hit_rank="{0}" peptide="{1}" peptide_prev_aa="" peptide_next_aa="" protein="{2}" num_tot_proteins="{3}" num_matched_ions="0" tot_num_ions="0" calc_neutral_pep_mass="{4}" massdiff="{5}" num_tol_term="2" num_missed_cleavages="{6}" num_matched_peptides="0">\n'.format(p.rank,p.OriginalPeptide,p.Proteins[0],p.ProteinCount,p.CalculatedParentMass,p.Massdiff,num_missing_cleavages)
                if p.ProteinCount>1:
                    for pro in p.Proteins[1:]:
                        psmItemString+='<alternative_protein protein="{0}"/>\n'.format(pro)
                if ('C' in p.IdentifiedPeptide) or ('~' in p.IdentifiedPeptide) or ('%' in p.IdentifiedPeptide):
                    if '%' not in p.IdentifiedPeptide:
                        psmItemString+='<modification_info modified_peptide="{0}">\n'.format(p.IdentifiedPeptide.replace('~','[147]'))
                    else:
                        psmItemString+='<modification_info modified_peptide="{0}" mod_nterm_mass="43.018425">\n'.format(p.IdentifiedPeptide.replace('~','[147]').replace('%','n[43]'))
                    for aa_id,aa in enumerate(p.IdentifiedPeptide):
                        if aa=='C':
                            psmItemString+='<mod_aminoacid_mass position="{0}" mass="160.030649" static="57.021464"/>\n'.format(str(aa_id+1-p.IdentifiedPeptide[:aa_id].count('~')))
                        if aa=='~':
                            psmItemString+='<mod_aminoacid_mass position="{0}" mass="147.035385" variable="15.994900" source="param"/>\n'.format(str(aa_id))
                    psmItemString+='</modification_info>\n'
                f.write(psmItemString)
                f.write('<search_score name="xcorr" value="0"/>\n<search_score name="deltacn" value="0"/>\n<search_score name="deltacnstar" value="0"/>\n<search_score name="spscore" value="0"/>\n<search_score name="sprank" value="0"/>\n<search_score name="expect" value="0"/>\n<analysis_result analysis="peptideprophet">\n')
                f.write('<peptideprophet_result probability="{0}" all_ntt_prob="(0.0000,0.0000,{0})">'.format(p.score))
                f.write('\n<search_score_summary>\n<parameter name="fval" value="0"/>\n<parameter name="ntt" value="2"/>\n<parameter name="nmc" value="0"/>\n<parameter name="massd" value="0.000"/>\n</search_score_summary>\n</peptideprophet_result>\n</analysis_result>\n</search_hit>\n')
            f.write('</search_result>\n')
            f.write('</spectrum_query>\n')
        f.write('</msms_run_summary>\n')
        f.write('</msms_pipeline_analysis>\n')
