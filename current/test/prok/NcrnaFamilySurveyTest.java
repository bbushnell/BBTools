package test.prok;

import java.io.File;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaFamilySurvey;
import prok.RrnaResourceIO;

/** Annotated locus identity, multipart territory and unequal-phylum denominator checks.
 * @author Brian Bushnell, Raiden */
public final class NcrnaFamilySurveyTest {
	public static void main(final String[] args){
		check(args.length==1, "Pass a fresh test directory");
		final File dir=new File(args[0]);check(dir.mkdir(), "Fixture directory must be new");
		check(new File(dir, "bacteria").mkdir() && new File(dir, "archaea").mkdir(), "Both domains must be represented");
		write(new File(dir, "bacteria/g1.gff.gz"), G1);
		write(new File(dir, "bacteria/g2.gff.gz"), "##gff-version 3\nc\tPGAP\tCDS\t1\t100\t.\t+\t0\tID=cds\n");
		write(new File(dir, "archaea/g3.gff.gz"), G3);
		write(new File(dir, "map.tsv"), "bacteria\tP1\tg1\nbacteria\tP1\tg2\narchaea\tP2\tg3\n");
		final long[] counters=new long[6];
		final ArrayList<NcrnaFamilySurvey.Locus> loci=NcrnaFamilySurvey.readLoci(new File(dir, "bacteria/g1.gff.gz").toString(), counters);
		check(loci.size()==14 && counters[4]==1, "Linked annotation layers and exact duplicate components collapse; neighboring copies survive");
		boolean labelKept=false, branchA=false, branchB=false, codingExcluded=false;
		boolean mixedRfam=false, mixedText=false, alpha=false, beta=false;
		int glycine=0;boolean multipart=false, overlapping=false;
		for(NcrnaFamilySurvey.Locus l:loci){
			if(l.family.equals("RF00504")){glycine++;multipart|=l.length==80 && l.spans.equals("1-40,61-100");}
			overlapping|=l.family.equals("RF00174") && l.length==51;
			labelKept|=l.family.equals("TEXT:16S ribosomal RNA") && l.length==50 && NcrnaFamilySurvey.shipped(l.family, l.type, l.label).equals("16S");
			branchA|=l.family.equals("RF99998") && l.spans.equals("1100-1119,1140-1159") && l.label.equals("branch A");
			branchB|=l.family.equals("RF99997") && l.spans.equals("1100-1104,1200-1209,1240-1249") && l.length==25 && l.label.equals("branch B");
			codingExcluded|=l.family.equals("RF99996") && l.length==10;
			mixedRfam|=l.family.equals("RF99995") && l.spans.equals("2000-2019,2040-2059") && l.label.equals("mixed Rfam");
			mixedText|=l.family.equals("TEXT:mixed text") && l.spans.equals("2100-2119,2140-2159") && l.length==40;
			alpha|=l.family.equals("TEXT:alpha") && l.spans.equals("2300-2319");
			beta|=l.family.equals("TEXT:beta") && l.spans.equals("2400-2419");
		}
		check(glycine==2 && multipart && overlapping, "Exon unions omit gaps and different overlapping families remain independent");
		check(labelKept && branchA && branchB && codingExcluded && counters[5]>0,
			"Unlabelled exons inherit family ancestry, preserve RNA labels, and cannot borrow a coding sibling's span");
		check(mixedRfam && mixedText && alpha && beta, "A shared gene parent cannot erase distinct Rfam/TEXT or TEXT/TEXT RNA branches");
		for(int threads:new int[]{1, 3}){
			final String prefix=new File(dir, "t"+threads).toString();
			NcrnaFamilySurvey.main(new String[]{"map="+new File(dir, "map.tsv"), "root="+dir, "expected=3", "outprefix="+prefix, "t="+threads});
			verify(prefix+".families.tsv");
			verifyAliases(prefix+".aliases.tsv");
			verifyExclusions(prefix+".excluded_regulatory.tsv");
			verifySummary(prefix+".summary.tsv");
		}
		compare(new File(dir, "t1.families.tsv"), new File(dir, "t3.families.tsv"));
		compare(new File(dir, "t1.phyla.tsv"), new File(dir, "t3.phyla.tsv"));
		compare(new File(dir, "t1.summary.tsv"), new File(dir, "t3.summary.tsv"));
		compare(new File(dir, "t1.excluded_regulatory.tsv"), new File(dir, "t3.excluded_regulatory.tsv"));
		write(new File(dir, "bad.gff"), "c\tPGAP\tncRNA\t20\t10\t.\t+\t.\tID=bad;Dbxref=RFAM:RF00504\n");
		boolean rejected=false;
		try{NcrnaFamilySurvey.readLoci(new File(dir, "bad.gff").toString(), new long[6]);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Invalid feature ranges cannot silently alter prevalence or territory");
		check(NcrnaFamilySurvey.shipped("RF00373", "RNase_P_RNA", "archaeal RNA").equals("RNaseP")
			&& NcrnaFamilySurvey.shipped("RF01685", "ncRNA", "6S").equals("6S")
			&& NcrnaFamilySurvey.shipped("RF00504", "riboswitch", "glycine").equals("NA"), "Shipped aliases cannot reappear as discovery candidates");
		check(NcrnaFamilySurvey.shipped("RF00177", "rRNA", "SSU_rRNA_bacteria").equals("16S")
			&& NcrnaFamilySurvey.shipped("RF01959", "rRNA", "SSU_rRNA_archaea").equals("16S")
			&& NcrnaFamilySurvey.shipped("RF02542", "rRNA", "SSU_rRNA_microsporidia").equals("18S"),
			"The shipped-alias table must classify microsporidian nuclear SSU as 18S, not 16S");
		for(String family:new String[]{"16S", "23S", "5S"}){
			final String label="ribosomal RNA-"+family;
			check(NcrnaFamilySurvey.shipped("TEXT:"+label, "rRNA", label).equals(family),
				"Alternate annotated ribosomal products must remain in the shipped-family block: "+label);
		}
		System.out.println("NCRNA_SURVEY_TEST_PASS layers multipart nearby overlap fallback zero_genome balanced median parallel malformed aliases regulatory_exclusions");
	}
	private static void verifyExclusions(final String path){
		final String[] classes={"minus_10_signal", "minus_35_signal", "other", "promoter", "ribosome_binding_site", "terminator"};
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "regulatory_class\tfeatures_excluded_without_Rfam_or_RNA_evidence", path);
			for(String type:classes){
				final byte[] line=in.nextLine();check(line!=null, "Each non-RNA regulatory class must remain visible in the exclusion ledger");p.set(line);
				check(p.terms()==2 && p.termEquals(type, 0) && p.parseInt(1)==(type.equals("minus_10_signal") ? 2 : 1),
					"Exclude one motif per class, plus a second genome's minus_10_signal; worker merging must retain both: "+type);
			}
			check(in.nextLine()==null, "A text-only riboswitch must survive rather than enter the non-RNA exclusion ledger");
		}finally{check(!in.close(), "Exclusion ledger must close without a read error");}
	}
	private static void verifySummary(final String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			check(in.nextLine()!=null, "Summary needs its header");final byte[] line=in.nextLine();check(line!=null, "Summary needs its measured row");p.set(line);
			check(p.terms()==11 && p.parseInt(2)==14 && p.parseInt(8)==7 && p.parseInt(9)==6 && p.parseInt(10)==1,
				"The fixture has7 Rfam identities,6 text keys and1 unknown; text labels are not biological families");
			check(in.nextLine()==null, "Summary contains exactly one corpus row");
		}finally{check(!in.close(), "Summary validation must close without a read error");}
	}

	/** Explicit reviewed identities protect against shifts in the paired alias arrays. */
	private static void verifyAliases(final String path){
		final String[] expected={"RF00177\t16S", "RF01959\t16S", "RF02540\t23S", "RF02541\t23S",
			"RF00001\t5S", "RF00005\ttRNA", "RF00010\tRNaseP", "RF00011\tRNaseP", "RF00373\tRNaseP",
			"RF00023\ttmRNA", "RF01849\ttmRNA", "RF01850\ttmRNA", "RF01851\ttmRNA", "RF04321\ttmRNA",
			"RF00013\t6S", "RF01685\t6S", "RF00169\tSRP", "RF04183\tSRP", "RF01854\tSRP", "RF01857\tSRP",
			"RF01852\ttRNA", "RF01960\t18S", "RF02542\t18S", "RF00002\t5.8S", "RF02543\teukLSU"};
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "rfam\tshipped_family", path);
			for(String row:expected){
				final byte[] line=in.nextLine();
				check(line!=null && row.equals(new String(line, java.nio.charset.StandardCharsets.UTF_8)),
					"Every emitted alias must match the independently reviewed identity: "+row);
				p.set(line);check(NcrnaFamilySurvey.shipped(p.parseString(0), "ncRNA", "unlabelled").equals(p.parseString(1)),
					"Discovery classification must agree with the emitted exclusion ledger: "+row);
			}
			check(in.nextLine()==null, "The audited exclusion ledger contains exactly25 Rfam identities");
		}finally{check(!in.close(), "Alias ledger validation must not ignore a read failure");}
	}
	private static void verify(final String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int rows=0;boolean found=false, shipped=false;
		try{
			check(in.nextLine()!=null, "Family output needs its schema");
			for(byte[] line; (line=in.nextLine())!=null;){
				p.set(line);check(p.terms()==16, "All rankable columns must be present");rows++;
				if(p.termEquals("RF00504", 1)){
					found=true;check(p.termEquals("CANDIDATE", 0) && p.parseLong(4)==110 && p.parseInt(6)==2 && p.parseInt(7)==3
						&& Math.abs(p.parseDouble(8)-2.0/3)<1e-8 && p.parseDouble(9)==.75 && p.parseInt(10)==2
						&& p.parseDouble(12)==20 && p.parseInt(13)==3 && p.parseDouble(14)==1,
						"RF00504 has lengths80/20/10 on2of3genomes: raw2/3, balanced(1/2+1)/2, median20");
				}
				if(p.termEquals("RF00373", 1)){shipped=p.termEquals("SHIPPED", 0) && p.termEquals("RNaseP", 3);}
				if(p.termEquals("TEXT:RNA switch", 1)){check(p.termEquals("UNRESOLVED", 0) && p.parseInt(13)==2,
					"Type-only and class-only riboswitches survive but their shared text label does not establish a Rfam identity");}
			}
		}finally{check(!in.close(), "Fixture result reader must close");}
		check(rows==14 && found && shipped, "All fourteen families, including unknown and shipped types, must survive");
	}
	private static void write(final File file, final String content){
		//TODO: Inherited fixture helper starts before protected ownership and closes
		//only after successful printing; annotate/report only in this survey pass.
		final ByteStreamWriter out=RrnaResourceIO.writer(file.toString());out.print(content);RrnaResourceIO.finish(out);
	}
	private static void compare(final File a, final File b){
		//TODO: Inherited paired-reader cleanup: second open can leak x; a throwing
		//first close can skip y or mask comparison failure. Both flags are checked
		//when closes return normally. Separate from production survey output repair.
		final ByteFile x=RrnaResourceIO.open(a.toString()), y=RrnaResourceIO.open(b.toString());
		try{for(byte[] line; (line=x.nextLine())!=null;){check(java.util.Arrays.equals(line, y.nextLine()), "Worker count must not change aggregate rows");}check(y.nextLine()==null, "Parallel output has extra rows");}
		finally{boolean error=x.close();error|=y.close();check(!error, "Both fixture comparison inputs must close");}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}

	private static final String G1="##gff-version 3\n"
		+"c\tPGAP\tgene\t1\t200\t.\t+\t.\tID=g;gene_biotype=ncRNA\n"
		+"c\tPGAP\tprimary_transcript\t1\t150\t.\t+\t.\tID=tx;Parent=g;gbkey=ncRNA;Dbxref=RFAM:RF00504;product=glycine\n"
		+"c\tPGAP\tncRNA\t1\t100\t.\t+\t.\tID=r;Parent=tx;Dbxref=RFAM:RF00504;product=glycine\n"
		+"c\tPGAP\texon\t1\t40\t.\t+\t.\tID=e1;Parent=r\n"
		+"c\tPGAP\texon\t61\t100\t.\t+\t.\tID=e2;Parent=r\n"
		+"c\tPGAP\tncRNA\t1\t100\t.\t+\t.\tID=dup;Dbxref=RFAM:RF00504;product=glycine\n"
		+"c\tPGAP\texon\t1\t40\t.\t+\t.\tParent=dup\n"
		+"c\tPGAP\texon\t61\t100\t.\t+\t.\tParent=dup\n"
		+"c\tPGAP\triboswitch\t220\t239\t.\t+\t.\tID=second;Dbxref=RFAM:RF00504;Note=glycine\n"
		+"c\tPGAP\tregulatory_region\t20\t70\t.\t+\t.\tID=overlap;Dbxref=RFAM:RF00174;Note=cobalamin\n"
		+"c\tPGAP\tregulatory_region\t300\t309\t.\t+\t.\tID=fallback;regulatory_class=riboswitch;Note=RNA%20switch\n"
		+"c\tPGAP\tminus_10_signal\t310\t318\t.\t+\t.\tID=dna10;gbkey=regulatory;regulatory_class=minus_10_signal;Note=putative promoter motif\n"
		+"c\tPGAP\tminus_35_signal\t320\t328\t.\t+\t.\tID=dna35;gbkey=regulatory;regulatory_class=minus_35_signal;Note=putative promoter motif\n"
		+"c\tPGAP\tregulatory_region\t330\t338\t.\t+\t.\tID=prom;gbkey=regulatory;regulatory_class=promoter;Note=promoter\n"
		+"c\tPGAP\tribosome_entry_site\t340\t345\t.\t+\t.\tID=rbs;gbkey=regulatory;regulatory_class=ribosome_binding_site;Note=putative RBS\n"
		+"c\tPGAP\tregulatory_region\t350\t359\t.\t+\t.\tID=term;gbkey=regulatory;regulatory_class=terminator;Note=terminator\n"
		+"c\tPGAP\triboswitch\t360\t369\t.\t+\t.\tID=typeOnly;Note=RNA%20switch\n"
		+"c\tPGAP\tregulatory_region\t370\t379\t.\t+\t.\tID=unclassified;gbkey=regulatory;regulatory_class=other;Note=unclassified regulatory feature\n"
		+"c\tPGAP\tncRNA\t400\t404\t.\t+\t.\tID=unknown\n"
		+"c\tPGAP\tgene\t500\t999\t.\t+\t.\tID=protein;gene_biotype=protein_coding\n"
		+"c\tPGAP\tCDS\t500\t999\t.\t+\t0\tID=cds;Parent=protein;product=enzyme\n"
		+"c\tPGAP\tmRNA\t500\t999\t.\t+\t.\tID=message;Parent=protein;product=enzyme\n"
		+"c\tPGAP\texon\t500\t999\t.\t+\t.\tParent=message;gbkey=CDS\n"
		+"c\tPGAP\tgene\t1000\t1049\t.\t+\t.\tID=legacyGene;gene_biotype=rRNA;product=16S ribosomal RNA\n"
		+"c\tPGAP\trRNA\t1000\t1049\t.\t+\t.\tID=legacy;Parent=legacyGene\n"
		+"c\tPGAP\texon\t1000\t1049\t.\t+\t.\tParent=legacy;gbkey=rRNA\n"
		+"c\tPGAP\tgene\t1100\t1299\t.\t+\t.\tID=multi;gene_biotype=ncRNA\n"
		+"c\tPGAP\tncRNA\t1100\t1159\t.\t+\t.\tID=branchA;Parent=multi;Dbxref=RFAM:RF99998;product=branch A\n"
		+"c\tPGAP\texon\t1100\t1119\t.\t+\t.\tParent=branchA\n"
		+"c\tPGAP\texon\t1140\t1159\t.\t+\t.\tParent=branchA\n"
		+"c\tPGAP\tncRNA\t1100\t1249\t.\t+\t.\tID=branchB;Parent=multi;Dbxref=RFAM:RF99997;product=branch B\n"
		+"c\tPGAP\texon\t1200\t1209\t.\t+\t.\tParent=branchB\n"
		+"c\tPGAP\texon\t1240\t1249\t.\t+\t.\tParent=branchB\n"
		+"c\tPGAP\texon\t1100\t1104\t.\t+\t.\tParent=branchA,branchB\n"
		+"c\tPGAP\tgene\t1400\t1999\t.\t+\t.\tID=coding;gene_biotype=protein_coding\n"
		+"c\tPGAP\tmRNA\t1400\t1999\t.\t+\t.\tID=codingRNA;Parent=coding\n"
		+"c\tPGAP\texon\t1400\t1999\t.\t+\t.\tParent=codingRNA\n"
		+"c\tPGAP\triboswitch\t1500\t1509\t.\t+\t.\tID=reg;Parent=coding;Dbxref=RFAM:RF99996;Note=coding neighbor\n"
		+"c\tPGAP\texon\t1500\t1504\t.\t+\t.\tParent=reg,codingRNA\n"
		+"c\tPGAP\texon\t1505\t1509\t.\t+\t.\tParent=reg\n"
		+"c\tPGAP\tgene\t2000\t2199\t.\t+\t.\tID=mixed;gene_biotype=ncRNA\n"
		+"c\tPGAP\tncRNA\t2000\t2059\t.\t+\t.\tID=mixedRfam;Parent=mixed;Dbxref=RFAM:RF99995;product=mixed Rfam\n"
		+"c\tPGAP\texon\t2000\t2019\t.\t+\t.\tParent=mixedRfam\n"
		+"c\tPGAP\texon\t2040\t2059\t.\t+\t.\tParent=mixedRfam\n"
		+"c\tPGAP\tncRNA\t2100\t2159\t.\t+\t.\tID=mixedText;Parent=mixed;product=mixed text\n"
		+"c\tPGAP\texon\t2100\t2119\t.\t+\t.\tParent=mixedText\n"
		+"c\tPGAP\texon\t2140\t2159\t.\t+\t.\tParent=mixedText\n"
		+"c\tPGAP\tgene\t2300\t2499\t.\t+\t.\tID=texts;gene_biotype=ncRNA\n"
		+"c\tPGAP\tncRNA\t2300\t2319\t.\t+\t.\tID=alpha;Parent=texts;product=alpha\n"
		+"c\tPGAP\tncRNA\t2400\t2419\t.\t+\t.\tID=beta;Parent=texts;product=beta\n";
	private static final String G3="##gff-version 3\n"
		+"a\tPGAP\triboswitch\t5\t14\t.\t-\t.\tID=reg;Dbxref=RFAM:RF00504;Note=glycine\n"
		+"a\tPGAP\tRNase_P_RNA\t20\t49\t.\t+\t.\tID=rnp;Dbxref=RFAM:RF00373;product=archaeal RNA\n"
		+"a\tPGAP\ttRNA\t60\t71\t.\t+\t.\tID=trna;product=tRNA-Leu\n"
		+"a\tPGAP\tminus_10_signal\t80\t88\t.\t+\t.\tID=secondDna;gbkey=regulatory;regulatory_class=minus_10_signal;Note=promoter motif\n";
}
