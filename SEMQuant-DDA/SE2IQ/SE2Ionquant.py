import sys
import re
import glob
db_type='metagenome-assembled' #uniprot or metagenome-assembled
fasta_path=sys.argv[1]
folder=sys.argv[2]
pro2psm_file=glob.glob(folder+'/*.pro2psm.txt')[0]
tab_file=glob.glob(folder+'/*.tab')[0]
extended_peptide_window=8
decoy_prefix='Rev2_'
psm_str='Spectrum\tSpectrum File\tPeptide\tModified Peptide\tExtended Peptide\tPrev AA\tNext AA\tPeptide Length\tCharge\tRetention\tObserved Mass\tCalibrated Observed Mass\tObserved M/Z\t' \
        'Calibrated Observed M/Z\tCalculated Peptide Mass\tCalculated M/Z\tDelta Mass\tExpectation\tHyperscore\tNextscore\tPeptideProphet Probability\tNumber of Enzymatic Termini\tNumber of Missed Cleavages' \
        '\tProtein Start\tProtein End\tIntensity\tAssigned Modifications\tObserved Modifications\tPurity\tIs ' \
        'Unique\tProtein\tProtein ID\tEntry Name\tGene\tProtein Description\tMapped Genes\tMapped Proteins\n'

pro_str='Protein\tProtein ID\tEntry Name\tGene\tLength\tOrganism\tProtein Description\tProtein ' \
        'Existence\tCoverage\tProtein Probability\tTop Peptide Probability\tTotal Peptides\tUnique Peptides\tRazor ' \
        'Peptides\tTotal Spectral Count\tUnique Spectral Count\tRazor Spectral Count\tTotal Intensity\tUnique ' \
        'Intensity\tRazor Intensity\tRazor Assigned Modifications\tRazor Observed Modifications\tIndistinguishable ' \
        'Proteins\n'

class PSM:
    def __init__(self, spectrum,filename, ParentCharge, MeasuredParentMass,CalculatedParentMass, Massdiff, score, peptide, modified_peptide, peptide_length, num_miss_cleavages, Retention, Protein,ProteinID, entry_name, gene, proteinDesc,mappedGenes, mappedProteins):
        self.spectrum=spectrum
        self.filename=filename
        self.ParentCharge = ParentCharge
        self.MeasuredParentMass=MeasuredParentMass
        self.CalculatedParentMass=CalculatedParentMass
        self.Massdiff = Massdiff
        self.score = score
        self.peptide=peptide
        self.modified_peptide = modified_peptide
        self.extended_peptide=''
        self.peptide_length=peptide_length
        self.retention=Retention
        self.num_miss_cleavages=num_miss_cleavages
        self.protein_start='0'
        self.protein_end='0'
        self.assigned_mod=[]
        self.is_unique=''
        self.Protein = Protein
        self.ProteinID = ProteinID
        self.entry_name = entry_name
        self.gene=gene
        self.proteinDesc=proteinDesc
        self.mappedGenes = mappedGenes
        self.mappedProteins=mappedProteins
    def update(self,seq):
        seq=seq.split('X')
        seq=''.join(seq)
        start_index=seq.index(self.peptide)
        self.protein_start=str(start_index+1)
        self.protein_end=str(start_index+int(self.peptide_length)+1)
        pre_pep=seq[start_index-8:start_index]
        post_pep=seq[start_index+int(self.peptide_length):start_index+int(self.peptide_length)+8]
        self.extended_peptide='.'.join([pre_pep,self.peptide,post_pep])
        temp_pep=self.modified_peptide
        if self.modified_peptide[0]=='%':
            self.assigned_mod.append('N-term(42.0106)')
            temp_pep=self.modified_peptide[1:]
        num_ptm=0
        for i in range(len(temp_pep)):
            if temp_pep[i]=='C':
                self.assigned_mod.append(str(i+1)+'C(57.0215)')
            elif temp_pep[i]=='~':
                self.assigned_mod.append(str(i-num_ptm)+'M(15.9949)')
                num_ptm+=1
        if len(self.mappedProteins)>0:
            self.is_unique='false'
        else:
            self.is_unique='true'


class Protein:
    def __init__(self,protein,proteinID,entry_name,gene,length,organism,proteinDesc,Indistinguishable_proteins):
        self.protein = protein
        self.proteinID = proteinID
        self.entry_name = entry_name
        self.gene = gene
        self.length = length
        self.organism = organism
        self.proteinDesc = proteinDesc
        self.coverage = 0.0
        self.protein_prob=0.0
        self.top_pep_prob=0.0
        self.modifications=[]
        self.indistinguishable_proteins=Indistinguishable_proteins
        self.PSMs=[]
    def update(self,seq):
        seqFlag=[0 for aa in seq]
        for psm in self.PSMs:
            for m in re.finditer(psm.peptide,seq):
                for i in range(m.start(),m.end()):
                    seqFlag[i]=1
        self.coverage=str(100*sum(seqFlag)/len(seq))
        pep_prob={}
        for psm in self.PSMs:
            if psm.modified_peptide not in pep_prob:
                pep_prob[psm.modified_peptide]=psm.score
            else:
                if psm.score>pep_prob[psm.modified_peptide]:
                    pep_prob[psm.modified_peptide]=psm.score
        protein_error_prob=1
        for pep in pep_prob:
            protein_error_prob*=(1-float(pep_prob[pep]))
        self.protein_prob=str(1-protein_error_prob)
        self.top_pep_prob=self.PSMs[0].score
        for psm in self.PSMs:
            self.modifications.extend(psm.assigned_mod)


def read_fasta(fasta_path):
    proID2seq={}
    with open(fasta_path) as f:
        for line in f:
            if line[0]=='>':
                id=line[1:].strip().split()[0]
            else:
                seq=line.strip()
                proID2seq[id]=seq
    return proID2seq

def read_tab(tab_file):
    scan2RT={}
    with open(tab_file) as f:
        for line_id, line in enumerate(f):
            if line_id>0:
                s=line.strip().split('\t')
                scan2RT[s[1]]=str(float(s[-2])*60) # transfer from mins to secs
    return scan2RT

def read_pro2psm(pro2psm_file,proID2seq,scan2RT):
    proteindict={}
    psmdict={}
    proteinDecoy = True
    with open(pro2psm_file) as f:
        for line_id, line in enumerate(f):
            if line_id>80:
                s=line.strip().split('\t')
                if s[0]=='+':
                    if s[-1]=='T':
                        proteinDecoy=False
                        protein_group=s[1].replace('{','').replace('}','').split(',')
                        target_protein_group=[]
                        for protein_item in protein_group:
                            if not protein_item.startswith(decoy_prefix):
                                target_protein_group.append(protein_item)
                        proteinDesc=s[8].replace('{','').replace('}','').split(',')
                        protein=target_protein_group[0]
                        for protein_item in target_protein_group:
                            if len(proID2seq[protein_item])>=len(proID2seq[protein]):
                                protein=protein_item
                        if db_type=='uniprot':
                            idx=protein.split('|')
                            proteinID=idx[1]
                            entry_name=idx[2]
                            gene=entry_name.split('_')[0]
                            organism=entry_name.split('_')[1]
                        else:
                            proteinID=''
                            entry_name=''
                            gene=''
                            organism=''
                        proteinLength=str(len(proID2seq[protein]))
                        indistinguishable_proteins = []
                        indistinguishable_genes =[]
                        if len(target_protein_group)>1:
                            index=target_protein_group.index(protein)
                            for i in range(len(target_protein_group)):
                                if i!=index:
                                    indistinguishable_proteins.append(target_protein_group[i])
                            proteinDesc=proteinDesc[index]
                            if db_type=='uniprot':
                                indistinguishable_genes=[x.strip().split('|')[1] for x in indistinguishable_proteins]
                            else:
                                indistinguishable_genes=[]
                        else:
                            proteinDesc=proteinDesc[0]
                        proteindict[protein]=Protein(protein,proteinID,entry_name,gene,proteinLength,organism,proteinDesc,indistinguishable_proteins)
                    else:
                        proteinDecoy = True
                elif s[0]=='*':
                    if proteinDecoy:
                        continue
                    else:
                        filename=s[1].split('.')[0]+'.pep.xml'
                        scan=s[2]
                        ParentCharge=s[3]
                        spectrum='.'.join([filename.split('.')[0],scan.zfill(5),scan.zfill(5),ParentCharge])
                        MeasuredParentMass=s[4]
                        CalculatedParentMass=s[5]
                        massdiff=str(float(MeasuredParentMass)-float(CalculatedParentMass))
                        score=s[11]
                        peptide=s[15][1:-1]
                        modified_peptide = s[14][1:-1]
                        peplength=str(len(s[15][1:-1]))
                        num_miss_cleavages=str(s[14][1:-2].count('K')+s[14][1:-2].count('R'))
                        Retention=scan2RT[scan]
                        if spectrum in psmdict:
                            psmdict[spectrum].mappedProteins.extend(indistinguishable_proteins)
                            psmdict[spectrum].mappedGenes.extend(indistinguishable_genes)
                            psmdict[spectrum].mappedProteins=list(set(psmdict[spectrum].mappedProteins))
                            psmdict[spectrum].mappedGenes=list(set(psmdict[spectrum].mappedGenes))
                        else:
                            psmdict[spectrum]=PSM(spectrum,filename, ParentCharge, MeasuredParentMass,CalculatedParentMass, massdiff, score, peptide, modified_peptide, peplength, num_miss_cleavages, Retention, protein,proteinID, entry_name, gene, proteinDesc,indistinguishable_genes, indistinguishable_proteins)
                        proteindict[protein].PSMs.append(psmdict[spectrum])

    return proteindict, psmdict

def writeout(proteindict, psmdict):
    with open(folder+'/psm.tsv','w') as f:
        f.write(psm_str)
        for key in psmdict:
            string=[key]
            string.append(psmdict[key].filename)
            string.append(psmdict[key].peptide)
            string.append(psmdict[key].modified_peptide.replace('~','[147]').replace('%','n[43]'))
            string.append(psmdict[key].extended_peptide)
            pre_pep=psmdict[key].extended_peptide.split('.')[0]
            post_pep=psmdict[key].extended_peptide.split('.')[2]
            if len(pre_pep)==0:
                string.append('-')
            else:
                string.append(pre_pep[-1])
            if len(post_pep)==0:
                string.append('-')
            else:
                string.append(post_pep[0])
            string.append(psmdict[key].peptide_length)
            string.append(psmdict[key].ParentCharge)
            string.append(psmdict[key].retention)
            string.append(psmdict[key].MeasuredParentMass)
            string.append(psmdict[key].MeasuredParentMass)
            string.append(str((float(psmdict[key].MeasuredParentMass)/int(psmdict[key].ParentCharge))+1.00728))
            string.append(str((float(psmdict[key].MeasuredParentMass)/int(psmdict[key].ParentCharge))+1.00728))
            string.append(psmdict[key].CalculatedParentMass)
            string.append(str((float(psmdict[key].CalculatedParentMass) / int(psmdict[key].ParentCharge)) + 1.00728))
            string.append(psmdict[key].Massdiff)
            string.append('0')
            string.append('0')
            string.append('0')
            string.append(psmdict[key].score)
            string.append('2')
            string.append(psmdict[key].num_miss_cleavages)
            string.append(psmdict[key].protein_start)
            string.append(psmdict[key].protein_end)
            string.append('0')
            string.append(','.join(psmdict[key].assigned_mod))
            string.append('')
            string.append('0')
            string.append(psmdict[key].is_unique)
            string.append(psmdict[key].Protein)
            string.append(psmdict[key].ProteinID)
            string.append(psmdict[key].entry_name)
            string.append(psmdict[key].gene)
            string.append(psmdict[key].proteinDesc)
            string.append(','.join(psmdict[key].mappedGenes))
            string.append(','.join(psmdict[key].mappedProteins))
            f.write('\t'.join(string)+'\n')

    with open(folder+'/protein.tsv','w') as f:
        f.write(pro_str)
        for key in proteindict:
            string=[key]
            string.append(proteindict[key].proteinID)
            string.append(proteindict[key].entry_name)
            string.append(proteindict[key].gene)
            string.append(proteindict[key].length)
            string.append(proteindict[key].organism)
            string.append(proteindict[key].proteinDesc)
            string.append('')
            string.append(str(proteindict[key].coverage))
            string.append(proteindict[key].protein_prob)
            string.append(proteindict[key].top_pep_prob)
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append('0')
            string.append(','.join(proteindict[key].modifications))
            string.append('')
            string.append(','.join(proteindict[key].indistinguishable_proteins))
            f.write('\t'.join(string)+'\n')









if __name__ == "__main__":
    proID2seq=read_fasta(fasta_path)
    scan2RT=read_tab(tab_file)
    proteindict, psmdict=read_pro2psm(pro2psm_file,proID2seq,scan2RT)
    for key in psmdict:
        psmdict[key].update(proID2seq[psmdict[key].Protein])
    for key in proteindict:
        proteindict[key].update(proID2seq[proteindict[key].protein])
    writeout(proteindict,psmdict)
    print('Done')