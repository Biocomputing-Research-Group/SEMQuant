import sys
import pandas as pd
threshold=0.15
class PSM:
    def __init__(self, filename, file, scan, ParentCharge, MeasuredParentMass,CalculatedParentMass,Massdiff, score,  IdentifiedPeptide,OriginalPeptide,
                 Proteins,ProteinCount,prev_aa,next_aa):
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
        self.prev_aa=prev_aa
        self.next_aa=next_aa

class Spectrum:
    def __init__(self,spectrumID,scan,precursor_neutral_mass,charge,PSM,retention):
        self.spectrumID=spectrumID
        self.scan=scan
        self.precursor_neutral_mass=precursor_neutral_mass
        self.charge=charge
        self.PSMs=[PSM]
        self.retention = retention
def extractExpMass(filename):
    scan2pep_dict = dict()

    with open(filename) as f:
        ifFirstScan = True
        for line_id, line in enumerate(f):
            if line_id > 189:
                if line[0] == '+':
                    if ifFirstScan:
                        scan = line.strip()
                        pep_list = []
                        ifFirstScan = False
                    else:
                        scan2pep_dict[scan] = pep_list
                        scan = line.strip()
                        pep_list = []
                elif line[0] == '*':
                    pep_list.append(line.strip())
    scan2pep_dict[scan] = pep_list

    with open(filename.replace('.target.','.decoy.')) as f:
        ifFirstScan = True
        for line_id, line in enumerate(f):
            if line_id > 189:
                if line[0] == '+':
                    if ifFirstScan:
                        scan = line.strip()
                        pep_list = []
                        ifFirstScan = False
                    else:
                        if scan in scan2pep_dict:
                            scan2pep_dict[scan].extend(pep_list)
                        else:
                            scan2pep_dict[scan] = pep_list
                        scan = line.strip()
                        pep_list = []
                elif line[0] == '*':
                    pep_list.append(line.strip())
    if scan in scan2pep_dict:
        scan2pep_dict[scan].extend(pep_list)
    else:
        scan2pep_dict[scan] = pep_list

    dMass = dict()
    for scan in scan2pep_dict:
        sn=scan.split('\t')[2]
        for pep in scan2pep_dict[scan]:
            x=pep.split('\t')
            dMass[sn+'_'+x[1]]=x[-1]
    del scan2pep_dict
    return dMass

def sipros_parse(filename,dMass):
    target_file=filename.replace('.pin','.target.tsv')
    decoy_file=filename.replace('.pin','.decoy.tsv')
    poutdf=pd.read_csv(target_file,sep='\t')
    poutdf=pd.concat([poutdf,pd.read_csv(decoy_file,sep='\t')],ignore_index=True)
    poutdf = poutdf.rename(columns={'PSMId': 'SpecId'})
    pindf=pd.read_csv(filename,sep='\t')
    df=pd.merge(pindf,poutdf,on='SpecId')
    df=df[df['posterior_error_prob']<threshold]
    df=df.reset_index(drop=True)
    df=df.sort_values('posterior_error_prob')
    del pindf
    del poutdf
    spectrum_d=dict()
    for idx in range(len(df)):
        SpecId=df['SpecId'][idx].strip().split('.')
        filename = '.'.join(SpecId[:1])
        file = SpecId[0]
        scan = df['ScanNr'][idx]
        charge = df['parentCharges'][idx]
        MeasuredParentMass = float(dMass[str(scan)+'_'+df['Peptide'][idx]])
        CalculatedParentMass = df['ExpMass'][idx]
        massdiff = str(abs(MeasuredParentMass - CalculatedParentMass))
        score = str(1-float(df['posterior_error_prob'][idx]))
        IdentifiedPeptide = df['Peptide'][idx][df['Peptide'][idx].index('[')+1:df['Peptide'][idx].index(']')]
        OriginalPeptide = IdentifiedPeptide.replace('~', '')
        Proteins = df['Proteins'][idx].replace('{', '')
        Proteins = Proteins.replace('}', '')
        Proteins = Proteins.split(',')
        rt=str(df['retentiontime'][idx]*60)
        prev_aa=''
        next_aa=''
        if df['Peptide'][idx][0]=='[':
            prev_aa='-'
        else:
            prev_aa=df['Peptide'][idx][0]
        if df['Peptide'][idx][-1]==']':
            next_aa='-'
        else:
            next_aa=df['Peptide'][idx][-1]
        idx = file + '.' + str(scan) + '.' + str(scan) + '.' + str(charge)+'_'+str(MeasuredParentMass)
        PSM_ins = PSM(filename, file, scan, charge, MeasuredParentMass, CalculatedParentMass, massdiff, score,
                      IdentifiedPeptide, OriginalPeptide, Proteins, len(Proteins),prev_aa,next_aa)
        if idx in spectrum_d:
            PSM_ins.rank = len(spectrum_d[idx].PSMs) + 1
            spectrum_d[idx].PSMs.append(PSM_ins)
        else:
            PSM_ins.rank = 1
            Spectrum_ins = Spectrum(idx, scan, MeasuredParentMass, charge, PSM_ins,rt)
            spectrum_d[idx] = Spectrum_ins

    return spectrum_d




if __name__ == "__main__":
    dMass = extractExpMass(sys.argv[1])
    spectrum_d = sipros_parse(sys.argv[2],dMass)
    with open('ProteinProphetHeader_DIANN.txt') as f:
        header=f.read()
    with open(sys.argv[3]+'/'+sys.argv[4]+'.pep.xml','w') as f:
        f.write(header)
        for spectrum in spectrum_d:
            spectrumItem=spectrum_d[spectrum]
            spectrumItemString='<spectrum_query spectrum="{0}" start_scan="{1}" end_scan="{1}" precursor_neutral_mass="{2}" assumed_charge="{3}" index="0" retention_time_sec="{4}">\n'.format('_'.join(spectrumItem.spectrumID.split('_')[:-1]),spectrumItem.scan,spectrumItem.precursor_neutral_mass,spectrumItem.charge,spectrumItem.retention)+'<search_result>\n'
            f.write(spectrumItemString)
            PSM_element=[spectrumItem.PSMs[0]]
            for p in PSM_element:
                tests=p.OriginalPeptide[0:-1]
                num_missing_cleavages=tests.count('K')+tests.count('R')
                psmItemString='<search_hit hit_rank="{0}" peptide="{1}" peptide_prev_aa="{2}" peptide_next_aa="{3}" protein="{4}" num_tot_proteins="{5}" num_matched_ions="0" tot_num_ions="0" calc_neutral_pep_mass="{6}" massdiff="{7}" num_tol_term="2" num_missed_cleavages="{8}" num_matched_peptides="0">\n'.format(p.rank,p.OriginalPeptide,p.prev_aa,p.next_aa,p.Proteins[0],p.ProteinCount,p.CalculatedParentMass,p.Massdiff,num_missing_cleavages)
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
                f.write('<peptideprophet_result probability="{0}" all_ntt_prob="({0},{0},{0})">'.format(p.score))
                f.write('\n<search_score_summary>\n<parameter name="fval" value="0"/>\n<parameter name="ntt" value="2"/>\n<parameter name="nmc" value="{0}"/>\n<parameter name="massd" value="0.000"/>\n</search_score_summary>\n</peptideprophet_result>\n</analysis_result>\n</search_hit>\n'.format(num_missing_cleavages))
            f.write('</search_result>\n')
            f.write('</spectrum_query>\n')
        f.write('</msms_run_summary>\n')
        f.write('</msms_pipeline_analysis>\n')
