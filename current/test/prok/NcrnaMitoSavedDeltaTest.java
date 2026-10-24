package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaMitoSavedDelta;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Exercises gross losses hidden by equal net TP and prior-contig ambiguity.
 * @author Raiden */
public final class NcrnaMitoSavedDeltaTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Explicit fresh fixture directory");
		final String d=args[0].substring(7);check(new File(d).mkdir(),"Fresh output directory required");
		final String h="arm\tpanel\tassembly\tcontig\tcompartment\tview\trole\tid\tlabel\tstart\tstop\tstrand\toutcome\tpartner\td5\td3\tgeometry\treference_kind\n";
		final String a=row("baseline","s1","TRUTH","r1",100,200,"TP","b1")+row("C","s1","TRUTH","r1",100,200,"FN","NA")
			+row("baseline","s1","TRUTH","r2",500,600,"FN","NA")+row("C","s1","TRUTH","r2",500,600,"TP","c1")
			+row("baseline","s1","CALL","b1",100,200,"TP","r1")+row("C","s1","CALL","c1",500,600,"TP","r2")
			+row("C","s1","CALL","fp1",100,210,"APPARENT_FP","NA")+row("C","s2","CALL","fp2",10,30,"APPARENT_FP","NA");
		write(d+"/audit.tsv",h+a);
		write(d+"/summary.tsv","arm\tpanel\tcompartment\tview\tTP\tapparentFP\tFN\texact5\texact3\tbothExact\tduplicateUnmatchedCalls\tunscoredCalls\tboundaryEligibleTP\tbothEndsEligibleTP\tprecision\trecall\tbothExactFraction\texact5Fraction\texact3Fraction\n"
			+"baseline\tP2\tNUCLEAR\tSSU\t1\t0\t1\t0\t0\t0\t0\t0\t1\t1\t1\t0.5\t0\t0\t0\nC\tP2\tNUCLEAR\tSSU\t1\t2\t1\t0\t0\t0\t0\t0\t1\t1\t0.333\t0.5\t0\t0\t0\n");
		final ByteBuilder old=new ByteBuilder();
		for(int i=0;i<32;i++){if(i>0){old.tab();}old.append(i==0?"case_id":i==27?"origin":"field"+i);}old.nl();
		for(int i=0;i<55;i++){
			for(int j=0;j<32;j++){if(j>0){old.tab();}
				if(j==0){old.append("old").append(i);}else if(j==3){old.append("a");}else if(j==4){old.append(i==0?"s1":"unrelated");}
				else if(j==5){old.append(105+i*1000);}else if(j==6){old.append(205+i*1000);}else if(j==7){old.append("+");}
				else if(j==18){old.append("0.91");}else if(j==27){old.append("NUMT");}else{old.append("NA");}
			}old.nl();
		}write(d+"/old.tsv",old.toString());
		final String dh="change\tseqid\tsource\ttype\tstart\tstop\tscore\tstrand\tphase\tattributes\n";
		write(d+"/p2.tsv",dh+"REMOVED\ts1\tBBTools\trRNA\t100\t200\t50\t+\t0\t18S,model:euk18S,family:euk18S\nADDED\ts1\tBBTools\trRNA\t100\t210\t60\t+\t0\t16S,model:mitochondria_coverage_test,family:16S\n");write(d+"/p3.tsv",dh);
		NcrnaMitoSavedDelta.main(new String[]{"audit="+d+"/audit.tsv","summary="+d+"/summary.tsv","old="+d+"/old.tsv","p2delta="+d+"/p2.tsv","p3delta="+d+"/p3.tsv","outprefix="+d+"/out"});
		check(read(d+"/out.summary.tsv").contains("P2\tNUCLEAR\tSSU\t2\t1\t1\t1\t1\t0\t2\t2\t0\t0\t1\t0\t1"),"Equal TP totals must preserve one gross loss and one gain; new55 relations conserve two calls");
		check(read(d+"/out.losses.tsv").contains("\tr1\t100\t200\t+\tb1\t100\t200\tfalse\t"),"Retain exact baseline call binding for the lost reference");
		check(read(d+"/out.neighbors.tsv").contains("MITO_16S") && read(d+"/out.old_matches.tsv").contains("SAME_STRAND_OVERLAP\t101\ttrue\tNUMT"),"Spatial overlap retains tier and old origin without inventing new identity");
		write(d+"/duplicate.tsv",h+a+row("baseline","s1","TRUTH","r1",100,200,"TP","b1"));
		boolean rejected=false;try{new NcrnaMitoSavedDelta().audit(d+"/duplicate.tsv");}catch(IllegalArgumentException e){rejected=true;}
		check(rejected,"Duplicate reference identity must fail loudly instead of inflating losses");
		System.out.println("MITO_SAVED_DELTA_TEST_PASS net_zero_gross_loss_gain old_overlap_new_contig original_partner duplicate_rejected");
	}
	private static String row(String arm,String seq,String role,String id,int start,int stop,String outcome,String partner){
		return arm+"\tP2\ta\t"+seq+"\tNUCLEAR\tSSU\t"+role+'\t'+id+"\tSSU\t"+start+'\t'+stop+"\t+\t"+outcome+'\t'+partner+"\tNA\tNA\tINTERNAL\tFROZEN_LOCUS_TRUTH_V2\n";
	}
	private static void write(String path,String text){final ByteStreamWriter out=RrnaResourceIO.writer(path);out.print(text);RrnaResourceIO.finish(out);}
	private static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] row;(row=in.nextLine())!=null;){b.append(row).nl();}}finally{check(!in.close(),"Fixture output must read completely");}return new String(b.toBytes(),StandardCharsets.UTF_8);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
