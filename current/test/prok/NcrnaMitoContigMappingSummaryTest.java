package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import prok.NcrnaMitoContigMappingSummary;
import prok.NcrnaMitoCallContexts;

/** Source-coordinate, secondary-union and missing-tail regression fixtures.
 * @author Brian, Ganyu */
public final class NcrnaMitoContigMappingSummaryTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1,"Supply a fresh fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		final Path meta=write(dir,"contigs.tsv","assembly\tcontig\tlength\tACGT\tGC\tGC_fraction_ACGT\texported_id\tsource_fasta\tsource_sha80\toriginal_header\n"
			+"A\tchrA\t12\t12\t6\t0.5\tMCONTIG_1\tNA\tNA\tfixture A\n"
			+"B\tchrB\t3\t3\t1\t0.333\tMCONTIG_2\tNA\tNA\tfixture B\n"
			+"C_POS\tC_POS\t5\t5\t2\t0.4\tNA\tNA\tNA\tpositive control\n");
		final String first=row("MCONTIG_1_0-4",0,"R1",1,"1S3=1I","AAAAA")
			+row("MCONTIG_1_0-4",272,"R1",20,"1S2=1X1S","*")
			+row("MCONTIG_1_0-4",272,"R2",30,"1=4S","*");
		final String second=row("MCONTIG_1_5-9",16,"R1",40,"1S2=1D1X1S","AAAAA");
		final String tail=row("MCONTIG_1_10-11",4,"*",0,"*","AA");
		final String other=row("MCONTIG_2_0-2",4,"*",0,"*","AAA");
		final String header="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:R1\tLN:100\n@SQ\tSN:R2\tLN:100\n";
		final Path sam=write(dir,"mapping.sam",header+first+second+tail+other);
		final String out=dir.resolve("result").toString();new NcrnaMitoContigMappingSummary().run(sam.toString(),meta.toString(),out,null,5,3);
		final List<String> lines=Files.readAllLines(Paths.get(out+".contigs.tsv"),StandardCharsets.UTF_8);
		check(lines.size()==3,"Both mapped and fully unmapped contigs must remain");final String[] a=lines.get(1).split("\t"),b=lines.get(2).split("\t");
		check(a[5].equals("3") && a[6].equals("2") && a[7].equals("6") && a[8].equals("7"),"Secondary overlap unions count each source base once, preserving unmapped tail denominator");
		check(a[11].equals("5") && a[12].equals("1") && a[13].equals("1") && a[14].equals("1") && a[15].equals("3") && Double.parseDouble(a[16])==.625,"Only primary CIGAR counts contribute to aggregate match/edit identity");
		check(a[17].equals("1") && a[18].equals("1") && a[19].equals("4"),"Ambiguity and reported-hit cap are visible without claiming exhaustive ties");
		check(b[5].equals("1") && b[6].equals("0") && b[8].equals("0") && b[16].equals("NA"),"A short entirely unmapped contig remains with an undefined identity");
		final String detail=NcrnaMitoCrossLibraryTest.read(out+".alignments.tsv.gz");
		check(detail.contains("\t7-8,6-6\t"),"Reverse-strand clip offsets must return to original forward source coordinates");
		final List<String> refs=Files.readAllLines(Paths.get(out+".references.tsv"),StandardCharsets.UTF_8);
		check(refs.get(1).startsWith("MCONTIG_1\tR1\t6\t") && refs.get(2).startsWith("MCONTIG_1\tR2\t1\t"),"Per-reference coverage unions secondary placements rather than summing them");
		final Path controlSam=write(dir,"with_control.sam",header+first+second+tail+other+row("C_POS_0-4",0,"R1",60,"5=","AAAAA"));
		final String controlOut=dir.resolve("with_control").toString();new NcrnaMitoContigMappingSummary().run(controlSam.toString(),meta.toString(),controlOut,null,"C_POS",5,3);
		final List<String> controlRows=Files.readAllLines(Paths.get(controlOut+".contigs.tsv"),StandardCharsets.UTF_8);
		check(controlRows.size()==4 && controlRows.get(1).startsWith("C_POS\t") && controlRows.get(1).split("\t")[20].equals("EXPLICIT_REFERENCE_CONTROL"),"Only explicitly named nonexported metadata rows may enter as positive controls");
		for(String cigar:new String[]{"2=2M1X","5M"}){
			final String label=cigar.equals("5M")?"all_M":"mixed_M";
			final Path msam=write(dir,label+".sam",header+other+row("C_POS_0-4",16,"R1",60,cigar,"AANNA"));
			final String mout=dir.resolve(label).toString();new NcrnaMitoContigMappingSummary().run(msam.toString(),meta.toString(),mout,"MCONTIG_2","C_POS",5,3);
			final String[] mr=Files.readAllLines(Paths.get(mout+".contigs.tsv"),StandardCharsets.UTF_8).get(1).split("\t");
			check(mr[16].equals("NA") && mr[24].equals("5") && mr[25].equals("5"),"M bases preserve possible aligned extent without a false point identity");
			if(label.equals("all_M")){check(mr[7].equals("0") && mr[21].equals("5") && Double.parseDouble(mr[22])==0 && Double.parseDouble(mr[23])==1,"M-only CIGAR supplies no known-base support and identity bounds zero to one");}
			else{check(mr[7].equals("3") && mr[21].equals("2") && Double.parseDouble(mr[22])==.4 && Double.parseDouble(mr[23])==.8,"Mixed M cannot become matches or be silently omitted from identity uncertainty");}
		}
		bad(dir,meta,header+first+second+other,"missing_tail","Missing tile");
		bad(dir,meta,header+first+second+tail+other+tail,"repeated_group","Repeated/noncontiguous");
		bad(dir,meta,header+row("MCONTIG_1_0-4",0,"R1",1,"1H5=","AAAAA")+second+tail+other,"hard_clip","Unsupported hard clip");
		bad(dir,meta,header+first+row("MCONTIG_1_5-9",16,"R1",40,"4=","AAAAA")+tail+other,"bad_cigar_length","every original tile base");
		contexts(dir,header);
		indelCalls(dir,header,false);indelCalls(dir,header,true);
		System.out.println("MITO_CONTIG_MAPPING_SUMMARY_TEST_PASS reverse_clips secondary_union primary_identity unmapped_tail cap missing_duplicate_invalid context_offsets call_identity hard_identity_filter origin_bins call_internal_deletions call_boundary_deletions insertion_crossings_both_strands");
	}
	private static void contexts(Path dir,String header)throws Exception{
		final Path contexts=write(dir,"contexts.tsv",NcrnaMitoCallContexts.CONTEXT_HEADER+"\n"
			+"X\tSOURCE_X\tA\tchrX\t1000\t1009\t10\t10\t300000\t300000\tCALL_NEIGHBORHOOD\tNA\tNA\tNA\tNA\tNA\n"
			+"Y\tSOURCE_Y\tA\tchrY\t0\t9\t10\t10\t10\t10\tWHOLE_CONTIG\tNA\tNA\tNA\tNA\tNA\n"
			+"Z\tSOURCE_Z\tA\tchrZ\t0\t4\t5\t5\t5\t5\tWHOLE_CONTIG\tNA\tNA\tNA\tNA\tNA\n");
		final Path calls=write(dir,"context_calls.tsv",NcrnaMitoCallContexts.CALL_HEADER+"\n"
			+"C1\tX\tSOURCE_X\tA\tchrX\t1003\t1007\t-\t2\t6\tfixture1\tmodel\n"
			+"C2\tY\tSOURCE_Y\tA\tchrY\t1\t10\t+\t0\t9\tfixture2\tmodel\n"
			+"C3\tY\tSOURCE_Y\tA\tchrY\t1\t5\t+\t0\t4\tfixture3\tmodel\n"
			+"C4\tZ\tSOURCE_Z\tA\tchrZ\t1\t5\t+\t0\t4\tfixture4\tmodel\n");
		final Path sam=write(dir,"context.sam",header+row("X_0-4",0,"R1",10,"5=","AAAAA")+row("X_5-9",16,"R1",30,"4=1X","AAAAA")
			+row("Y_0-4",0,"R1",10,"3=2X","AAAAA")+row("Y_5-9",0,"R1",30,"4=1M","AAAAN")+row("Z_0-4",0,"R1",10,"5M","NNNNN"));
		final String out=dir.resolve("context_result").toString();new NcrnaMitoContigMappingSummary().run(sam.toString(),null,out,null,null,5,3,contexts.toString(),calls.toString());
		final List<String> lines=Files.readAllLines(Paths.get(out+".calls.tsv"),StandardCharsets.UTF_8);check(lines.size()==5,"Every call must remain, including below-threshold and all-M mappings");
		final String[] x=lines.get(1).split("\t"),y=lines.get(2).split("\t"),z=lines.get(3).split("\t"),m=lines.get(4).split("\t");
		check(x[5].equals("1003") && x[6].equals("1007") && x[10].equals("5") && x[13].equals("4") && x[14].equals("1") && Double.parseDouble(x[18])==.8 && x[27].equals("NUMT"),"Call identity must restrict each primary CIGAR to the call and project reverse tiles into original coordinates");
		check(y[10].equals("4") && y[11].equals("5") && y[17].equals("1") && y[18].equals("NA") && Double.parseDouble(y[19])==.8 && y[27].equals("AMBIGUOUS") && y[30].equals("5"),"Only passing whole alignments enter coverage; ambiguous M stays outside known coverage and subthreshold support is separate");
		check(z[10].equals("0") && z[27].equals("UNMAPPED") && z[30].equals("5") && m[10].equals("0") && m[27].equals("UNMAPPED"),"A perfect call subspan cannot override a failed whole-alignment gate; all-M lower bound is zero");
		final String align=NcrnaMitoCrossLibraryTest.read(out+".call_alignments.tsv.gz");check(align.contains("1006-1006,1005-1005\t\t33-33,34-34"),"Restricted reference coordinates must follow the reverse-strand CIGAR, not a proportional projection");
		check(NcrnaMitoContigMappingSummary.passesIdentity(76,24,0,0,0) && !NcrnaMitoContigMappingSummary.passesIdentity(75,24,0,0,1),"Hard measured .76 identity is inclusive and unknown M contributes only to the denominator");
		check(NcrnaMitoContigMappingSummary.classify(199999,true,.8,.8).equals("MITO_CONTIG") && NcrnaMitoContigMappingSummary.classify(200000,false,1,.8).equals("NUMT")
			&& NcrnaMitoContigMappingSummary.classify(100,true,.79,.8).equals("NUMT") && NcrnaMitoContigMappingSummary.classify(100,true,1,.2).equals("AMBIGUOUS")
			&& NcrnaMitoContigMappingSummary.classify(100,true,1,.19).equals("UNMAPPED"),"Descriptive origin boundaries must match the authorized 200kb and 80/20-percent rules exactly");
	}
	private static void bad(Path dir,Path meta,String sam,String label,String expected)throws Exception{
		final Path path=write(dir,label+".sam",sam);boolean failed=false;
		try{new NcrnaMitoContigMappingSummary().run(path.toString(),meta.toString(),dir.resolve(label).toString(),null,5,3);}
		catch(IllegalArgumentException e){failed=e.getMessage().contains(expected);}
		check(failed,"Invalid input must fail loudly: "+label);
	}
	private static void indelCalls(Path dir,String header,boolean reverse)throws Exception{
		final String label=reverse?"indels_reverse":"indels_forward";
		final Path contexts=write(dir,label+".contexts.tsv",NcrnaMitoCallContexts.CONTEXT_HEADER+"\n"
			+"D\tSOURCE_D\tA\tchrD\t100\t119\t20\t20\t300000\t300000\tCALL_NEIGHBORHOOD\tNA\tNA\tNA\tNA\tNA\n"
			+"I\tSOURCE_I\tA\tchrI\t100\t119\t20\t20\t300000\t300000\tCALL_NEIGHBORHOOD\tNA\tNA\tNA\tNA\tNA\n");
		// Query positions 0..4, D(2), 5..14, D(1), 15..19.
		// Calls exercise both interior gaps, both exterior gaps, then each boundary separately.
		final Path calls=write(dir,label+".calls.tsv",NcrnaMitoCallContexts.CALL_HEADER+"\n"
			+indelCall("D_inside","D",4,15,reverse)+indelCall("D_outside","D",5,14,reverse)
			+indelCall("D_left_inside","D",4,14,reverse)+indelCall("D_right_inside","D",5,15,reverse)
			// Insertion occupies q5..7. One call starts in it; another ends in it.
			+indelCall("I_start_inside","I",6,9,reverse)+indelCall("I_stop_inside","I",3,6,reverse)
			+indelCall("I_whole_gap","I",4,8,reverse));
		final String sequence="AAAAAAAAAAAAAAAAAAAA";final int flag=reverse?16:0;
		final Path sam=write(dir,label+".sam",header+row("D_0-19",flag,"R1",10,"5=2D10=1D5=",sequence)+row("I_0-19",flag,"R1",10,"5=3I12=",sequence));
		final String out=dir.resolve(label+"_result").toString();new NcrnaMitoContigMappingSummary().run(sam.toString(),null,out,null,null,20,3,contexts.toString(),calls.toString());
		final List<String> result=Files.readAllLines(Paths.get(out+".calls.tsv"),StandardCharsets.UTF_8);check(result.size()==8,"Every boundary fixture must survive the whole-alignment .76 gate");
		checkIndel(result.get(1),"D_inside",12,0,3,12,.8);
		checkIndel(result.get(2),"D_outside",10,0,0,10,1);
		checkIndel(result.get(3),"D_left_inside",11,0,2,11,11.0/13);
		checkIndel(result.get(4),"D_right_inside",11,0,1,11,11.0/12);
		checkIndel(result.get(5),"I_start_inside",2,2,0,2,.5);
		checkIndel(result.get(6),"I_stop_inside",2,2,0,2,.5);
		checkIndel(result.get(7),"I_whole_gap",2,3,0,2,.4);
		final String detail=NcrnaMitoCrossLibraryTest.read(out+".call_alignments.tsv.gz");
		check(detail.contains("\t14-14,17-26,28-28\t"),"Internal deleted reference bases must not become supported intervals on either strand");
		final String[] inserted=result.get(5).split("\t");
		check(inserted[5].equals(reverse?"111":"107") && inserted[6].equals(reverse?"114":"110"),"Insertion-crossing calls retain original forward-source coordinates on both strands");
	}
	private static String indelCall(String id,String context,int qlo,int qhi,boolean reverse){
		final int start=reverse?19-qhi:qlo,stop=reverse?19-qlo:qhi;
		return id+'\t'+context+"\tSOURCE_"+context+"\tA\tchr"+context+'\t'+(101+start)+'\t'+(101+stop)+'\t'+(reverse?"-":"+")+'\t'+start+'\t'+stop+"\tfixture\tmodel\n";
	}
	private static void checkIndel(String line,String id,int matches,int insertions,int deletions,int covered,double identity){
		final String[] p=line.split("\t");
		check(p[0].equals(id) && Integer.parseInt(p[13])==matches && Integer.parseInt(p[15])==insertions && Integer.parseInt(p[16])==deletions
			&& Integer.parseInt(p[9])==covered && Integer.parseInt(p[10])==covered && Math.abs(Double.parseDouble(p[18])-identity)<1e-9,
			"Call-restricted edit counts/coverage/identity disagree with the explicit query-position fixture: "+id+" row="+line);
	}
	private static String row(String name,int flag,String ref,int pos,String cigar,String seq){
		return name+'\t'+flag+'\t'+ref+'\t'+pos+"\t60\t"+cigar+"\t*\t0\t0\t"+seq+"\t*\tYI:f:75.00\tYR:i:100\n";
	}
	private static Path write(Path dir,String name,String text)throws Exception{final Path path=dir.resolve(name);Files.write(path,text.getBytes(StandardCharsets.UTF_8));return path;}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
