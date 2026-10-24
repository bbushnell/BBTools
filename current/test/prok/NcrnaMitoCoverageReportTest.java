package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import prok.NcrnaMitoCoverageReport;

/** Copy/unique conservation, cross-lineage aliases, dominant ties and bad evidence.
 * @author Brian, Ganyu */
public final class NcrnaMitoCoverageReportTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be fresh");
		write(root,"runs.tsv","cut\tprefix\n0.60\t"+root+"/low\n0.66\t"+root+"/high\n");
		final ArrayList<String> files=new ArrayList<String>();files.add("runs.tsv");
		for(String name:new String[]{"low","high"}){final boolean high=name.equals("high");final String selected=high?"nohit:B":"ssu:A";final String cut=high?"0.66":"0.60";
			write(root,name+".models.tsv",NcrnaMitoCoverageReport.MODEL_HEADER+"\nssu:A\tmitochondria_coverage_ssu_A\t1000\tssu:A\tssu\tA n=3\nnohit:B\tmitochondria_coverage_nohit_B\t800\tnohit:B\tnohit\tB n=2\n");
			write(root,name+".summary.tsv",summary(cut,high,0));
			write(root,name+".selection.tsv",NcrnaMitoCoverageReport.STEP_HEADER+"\n1\t"+selected+(high?"\t1\t2\t1\t2\n":"\t2\t3\t2\t1\n"));
			write(root,name+".members.tsv.gz",NcrnaMitoCoverageReport.MEMBER_HEADER+"\n"
				+member("CONSTRUCTION","X1","U0","A",false,selected,selected,high?.9f:.8f,"nohit:B",.9f,true)
				+member("CONSTRUCTION","X2","U0","B",false,selected,selected,high?.9f:.8f,"nohit:B",.9f,true)
				+member("CONSTRUCTION","X3","U1","A",false,high?"NA":selected,selected,high?.4f:.62f,"ssu:A",.62f,!high)
				+member("PANEL_INCLUDED","P1","U2","B",true,"NA",selected,high?.3f:.5f,"ssu:A",.5f,false));
			write(root,name+".strata.tsv",strata(high,false));
			write(root,name+".composition.tsv",NcrnaMitoCoverageReport.COMPOSITION_HEADER+"\n"
				+selected+"\tSOURCE_CLUSTER_FOR_SELECTED_"+selected+"\tA\t"+(high?"1\t1":"2\t2")+"\n"
				+selected+"\tSOURCE_CLUSTER_FOR_SELECTED_"+selected+"\tB\t1\t1\n"
				+selected+"\tNEWLY_COVERED_TARGETS\tA\t"+(high?"1\t1":"2\t2")+"\n"
				+selected+"\tNEWLY_COVERED_TARGETS\tB\t1\t1\n");
			for(String suffix:new String[]{".models.tsv",".summary.tsv",".selection.tsv",".members.tsv.gz",".strata.tsv",".composition.tsv"}){files.add(name+suffix);}
		}
		String pins=NcrnaMitoCrossLibraryTest.bind(root,sha,files.toArray(new String[0]));run(root,pins,"good");
		final String curve=read(root+"/good.curve.tsv");final String[] low=NcrnaMitoCrossLibraryTest.find(curve,"0.6\tALL\t"),high=NcrnaMitoCrossLibraryTest.find(curve,"0.66\tALL\t");
		check(low[3].equals("4") && low[4].equals("3") && low[5].equals("3") && low[6].equals("2") && high[5].equals("2") && high[6].equals("1") && high[11].equals("0"),"Threshold curve must preserve copies, unique sequences and zero new alignments");
		final String residual=read(root+"/good.residual_lineages.tsv");check(residual.contains("0.66\tALL\tALL\tALL\tALL\t4\t3\t2\t1\t2\t2\tfalse"),"Global unique count cannot be the sum of overlapping lineage memberships");
		check(residual.contains("0.66\tPHYLUM\tALL\tALL\tA\t2\t2\t1\t1\t1\t1") && residual.contains("0.66\tPHYLUM\tALL\tALL\tB\t2\t2\t1\t1\t1\t1"),"Residual lineage rows retain cross-lineage sequence copies without inventing extra global uniques");
		final String models=read(root+"/good.models.tsv");final String[] source=NcrnaMitoCrossLibraryTest.find(models,"0.66\t1\tnohit:B\tnohit:B\tmitochondria_coverage_nohit_B\t800\tSOURCE_CLUSTER\t"),targets=NcrnaMitoCrossLibraryTest.find(models,"0.66\t1\tnohit:B\tnohit:B\tmitochondria_coverage_nohit_B\t800\tNEWLY_COVERED_TARGETS\t");
		check(source[8].equals("NA") && source[9].equals("2") && targets[8].equals("1") && targets[9].equals("2"),"Source aggregate unique incidences cannot masquerade as global unique counts; assigned targets permit exact union");
		check(source[10].equals("A;B") && source[11].equals("2") && source[13].equals("A;B") && source[14].equals("2"),"Every dominant-lineage tie must be retained in deterministic order");
		write(root,"high.summary.tsv",summary("0.66",true,1));pins=NcrnaMitoCrossLibraryTest.bind(root,sha,files.toArray(new String[0]));bad(root,pins,"nonzero_alignments");
		write(root,"high.summary.tsv",summary("0.66",true,0));write(root,"high.strata.tsv",strata(true,true));pins=NcrnaMitoCrossLibraryTest.bind(root,sha,files.toArray(new String[0]));bad(root,pins,"inconsistent_stratum");
		System.out.println("MITO_COVERAGE_REPORT_TEST_PASS curve copy_unique residual_lineages overlapping_lineages source_vs_assigned explicit_ties zero_alignment_guard conservation_failure");
	}
	private static String summary(String cut,boolean high,int alignments){return NcrnaMitoCoverageReport.SUMMARY_HEADER+"\n4\t3\t2\t1\t"+cut+(high?"\t1\t2\t2\t2":"\t2\t1\t3\t1")+"\t1\t1\tSHORTER_QUERY_COVERAGE_NOT_FULL_MODEL_CALLER_VERIFICATION\tselect\t2\t"+alignments+"\t6\n";}
	private static String strata(boolean high,boolean bad){return NcrnaMitoCoverageReport.STRATA_HEADER+"\nCONSTRUCTION:NOT_MARKED_PARTIAL\t3\t2"+(high && !bad?"\t2\t1\t1\t1":"\t3\t2\t0\t0")+"\nPANEL_INCLUDED:PARTIAL\t1\t1\t0\t0\t1\t1\n";}
	private static String member(String role,String id,String uid,String phylum,boolean partial,String assigned,String best,float score,String pool,float poolScore,boolean covered){final String[] p=new String[19];Arrays.fill(p,"NA");p[0]=role;p[1]=id;p[2]=uid;p[3]=phylum.equals("A")?"4":"5";p[4]=phylum;p[5]=""+partial;p[6]=role.equals("PANEL_INCLUDED")?"NA":"false";p[7]=partial?"PARTIAL_INTERNAL":"ssu";p[9]=assigned;p[10]=best;p[11]=""+score;p[12]=""+score;p[13]="1";p[14]="1";p[15]=pool;p[16]=""+poolScore;p[17]=covered?"COVERED_SHORTER_QUERY":"UNCOVERED_AT_CUT";p[18]="fixture_source";return String.join("\t",p)+"\n";}
	private static void run(String root,String pins,String out){NcrnaMitoCoverageReport.main(new String[]{"runs="+root+"/runs.tsv","verified="+pins,"outprefix="+root+"/"+out,"cuts=.60,.66","target=.66"});}
	private static void bad(String root,String pins,String label){boolean failed=false;try{run(root,pins,label);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+"/"+label+".curve.tsv").exists(),"Inconsistent evidence must fail before report emission: "+label);}
	private static void write(String root,String name,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+name,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
