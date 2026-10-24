package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.EukPanelTruthFreeze;
import prok.NcrnaClippingComparison;
import prok.NcrnaClippingComparison.Genome;
import prok.NcrnaClippingComparison.Contig;
import prok.NcrnaClippingComparison.Feature;
import structures.ByteBuilder;

/** Frozen truth changes reference units, not the one-to-one or duplicate rules.
 * @author Raiden */
public final class NcrnaFrozenLocusComparisonTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Fresh outdir required");final String root=args[0].substring(7);require(new File(root).mkdir(),"Fresh fixture root required");
		final String truth=EukPanelTruthFreeze.HEADER+'\n'+truth(1,100,199,true,true,"FULLSPAN_INTACT")+truth(2,300,399,false,false,"FULLSPAN_INTACT")+truth(3,600,699,true,false,"PARTIAL_INTERNAL");
		write(root+"/truth.tsv",truth);
		final String history=NcrnaClippingComparison.HISTORY_HEADER+'\n'+history(1,100,199,"+","LSU",true)+history(2,400,499,"+","LSU",true)
			+history(3,500,599,"+","5.8S",true)+history(4,700,799,"-","LSU",true)+history(5,800,899,"+","UNKNOWN",false)
			+history(6,800,899,"+","LSU",false).replace("\t800\t899\t","\tx\ty\t");
		write(root+"/history.tsv",history);
		final NcrnaClippingComparison x=new NcrnaClippingComparison();final Genome g=new Genome("A","P2","/unused.fa","00000000000000000000");x.addGenome(g);
		final Contig c=new Contig(g,"chr",1000,"NUCLEAR");g.contigs.put("chr",c);g.ordered.add(c);x.modelFamilies.put("RF02543","LSU");
		x.loadFrozenTruth(root+"/truth.tsv",3);x.loadHistoricalMasks(root+"/history.tsv",6);
		final int[] starts={100,100,300,400,500,700,600,800};for(int i=0;i<starts.length;i++){c.calls.add(new Feature(c,"LSU","lsu",starts[i],starts[i]+99,"+","Q"+i,"NA","CALL"));}
		final ByteStreamWriter audit=new ByteStreamWriter(root+"/audit.tsv",false,false,true);audit.start();x.grade(c,"TEST","LSU",audit,new ByteBuilder());require(!audit.poisonAndWait(),"Audit fixture write failed");
		final long[] actual=x.totals.get("TEST\tP2\tNUCLEAR\tLSU"),expected={2,4,0,2,2,1,1,2,2,1};
		require(Arrays.equals(actual,expected),"Expected TP/FP/FN/exact5/exact3/both/duplicate/unscored/boundaryTP/bothTP="+Arrays.toString(expected)+" got "+Arrays.toString(actual));
		final long[] intact=x.truthStrata.get("TEST\tP2\tNUCLEAR\tLSU\tFULLSPAN_INTACT");
		require(Arrays.equals(intact,new long[]{1,0,1,1,0,1,0}),"Crossfold truth remains detected and exact at each end; only both-end denominator is excluded");
		final String text=read(root+"/audit.tsv");require(text.contains("DUPLICATE_APPARENT_FP") && text.contains("UNSCORED_HIST_ONLY") && text.contains("UNSCORED_OVERLAP") && text.contains("FROZEN_LOCUS_TRUTH_V1"),"Audit must retain distinct duplicate, historical and model-mask outcomes and frozen provenance");
		write(root+"/unfrozen.tsv",truth.replace("\tTRUTH_V1\t","\tFREEZE_CANDIDATE_NOT_FROZEN\t"));
		boolean failed=false;try{new NcrnaClippingComparison().loadFrozenTruth(root+"/unfrozen.tsv",3);}catch(IllegalArgumentException e){require(e.getMessage().contains("selected frozen truth"),"Wrong rejection: "+e);failed=true;}require(failed,"An unfrozen candidate cannot silently become scored truth");
		final String v2=truth.replace("TRUTH_V1","TRUTH_V2").replace("1.5_WITH_HOMOLOG_RULE_WITHDRAWN",EukPanelTruthFreeze.V2_POLICY);write(root+"/v2.tsv",v2);
		final NcrnaClippingComparison x2=new NcrnaClippingComparison();final Genome g2=new Genome("A","P2","/unused.fa","00000000000000000000");x2.addGenome(g2);
		final Contig c2=new Contig(g2,"chr",1000,"NUCLEAR");g2.contigs.put("chr",c2);g2.ordered.add(c2);x2.modelFamilies.put("RF02543","LSU");
		x2.loadFrozenTruth(root+"/v2.tsv",3,EukPanelTruthFreeze.V2_POLICY,2);x2.loadHistoricalMasks(root+"/history.tsv",6);
		for(int i=0;i<starts.length;i++){c2.calls.add(new Feature(c2,"LSU","lsu",starts[i],starts[i]+99,"+","Q"+i,"NA","CALL"));}
		final ByteStreamWriter audit2=new ByteStreamWriter(root+"/audit2.tsv",false,false,true);audit2.start();x2.grade(c2,"TEST","LSU",audit2,new ByteBuilder());require(!audit2.poisonAndWait(),"Version2 audit write failed");
		require(Arrays.equals(actual,x2.totals.get("TEST\tP2\tNUCLEAR\tLSU")),"Version labeling must not change matching or boundary denominators");
		require(read(root+"/audit2.tsv").replace("FROZEN_LOCUS_TRUTH_V2","FROZEN_LOCUS_TRUTH_V1").equals(text),"Version2 audit must differ only in its explicit reference label");
		failed=false;try{new NcrnaClippingComparison().loadFrozenTruth(root+"/v2.tsv",3);}catch(IllegalArgumentException e){failed=e.getMessage().contains("selected frozen truth");}require(failed,"Version2 may not be silently loaded as version1");
		final String summaryHeader="arm\tpanel\tcompartment\tview\tTP\tapparentFP\tFN\texact5\texact3\tbothExact\tduplicateUnmatchedCalls\tunscoredCalls\tboundaryEligibleTP\tbothEndsEligibleTP\tprecision\trecall\tbothExactFraction\texact5Fraction\texact3Fraction\n";
		final String baselineRow="TEST\tP2\tNUCLEAR\tLSU\t2\t4\t0\t2\t2\t1\t1\t2\t2\t1\t0.33333333\t1\t1\t1\t1\n";
		write(root+"/baseline_summary.tsv",summaryHeader+baselineRow);
		x2.writeReferenceDelta(root+"/baseline_summary.tsv",root+"/equal");
		final long[] changed=x2.totals.get("TEST\tP2\tNUCLEAR\tLSU");changed[1]++;changed[7]--;
		x2.writeReferenceDelta(root+"/baseline_summary.tsv",root+"/changed");
		final ByteFile delta=ByteFile.makeByteFile(FileFormat.testInput(root+"/changed.reference_delta.tsv",FileFormat.TEXT,null,true,true));final parse.LineParser1 dp=new parse.LineParser1('\t');
		try{delta.nextLine();dp.set(delta.nextLine());require(dp.terms()==21 && dp.parseLong(5)==4 && dp.parseLong(8)==5 && dp.parseLong(11)==1 && dp.parseLong(13)==2 && dp.parseLong(14)==1,"Lifted mask must increase FP and reduce unscored");
			require(Math.abs(dp.parseDouble(15)-1.0/3)<1e-8 && Math.abs(dp.parseDouble(16)-2.0/7)<1e-8 && dp.parseDouble(17)==1 && dp.parseDouble(18)==1 && dp.parseDouble(19)==1 && dp.parseDouble(20)==1,"Precision changes while recall and exact-boundary fractions remain fixed");
		}finally{require(!delta.close(),"Delta fixture read failed");}
		x2.totals.clear();x2.totals.put("TEST\tP1\tPROKARYOTIC\tLSU",changed);
		write(root+"/prok_summary.tsv",summaryHeader+baselineRow.replace("P2\tNUCLEAR","P1\tPROKARYOTIC"));
		failed=false;try{x2.writeReferenceDelta(root+"/prok_summary.tsv",root+"/bad_prok_delta");}catch(IllegalArgumentException e){failed=e.getMessage().contains("untouched prokaryotic panel");}require(failed,"Euk reference changes must not silently change prok metrics");
		System.out.println("PASS frozen locus matching, both-only eligibility, below-GA/HIST masks, duplicate priority, wrong-family/strand and unclassified-source controls");
	}
	static String truth(int id,int start,int end,boolean admitted,boolean cross,String geometry){
		final String[] p=new String[35];Arrays.fill(p,"0");p[0]=""+id;p[1]="A";p[2]="chr";p[3]="+";p[4]=""+start;p[5]=""+end;p[6]="LSU";p[7]="test";p[8]="NUCLEAR";p[9]="RF02543";p[10]="run";p[11]="/raw.tbl";p[12]=""+id;p[13]="80";p[14]=admitted?"1":"0";p[15]=geometry.equals("FULLSPAN_INTACT")?"1":"0";p[17]="INTERNAL";p[18]="1";p[19]="RF02543";p[22]="NONE";p[23]="DESCRIPTIVE_ONLY";p[24]=admitted?"1":"0";p[25]=admitted?"RF02543":"NA";p[26]=admitted?"ADMITTED":"BELOW_GA";p[27]=geometry;p[28]=cross?"AMBIGUOUS_CROSSFOLD":"NONE";p[31]=admitted?"1":"0";p[32]="TRUTH_V1";p[33]=admitted && !cross?"1":"0";p[34]="1.5_WITH_HOMOLOG_RULE_WITHDRAWN";return String.join("\t",p)+'\n';
	}
	static String history(int id,int start,int end,String strand,String family,boolean eligible){return "A\tchr\t"+strand+'\t'+start+'\t'+end+'\t'+family+"\tHIST_ONLY\t"+(eligible?1:0)+"\t/source.gff\t00000000000000000000\t"+id+"\trRNA\tproduct=large subunit ribosomal RNA\t"+(eligible?"VALID_UNRESOLVED_SOURCE":"UNKNOWN_FAMILY")+'\n';}
	static void write(String path,String text){final ByteStreamWriter w=new ByteStreamWriter(path,false,false,true);w.start();w.print(text);require(!w.poisonAndWait(),"Fixture write failed");}
	static String read(String path){final ByteFile r=ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=r.nextLine())!=null;){b.append(line).nl();}}finally{require(!r.close(),"Fixture read failed");}return new String(b.toBytes(),StandardCharsets.UTF_8);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
