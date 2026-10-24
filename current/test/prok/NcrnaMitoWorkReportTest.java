package test.prok;

import java.io.File;
import prok.NcrnaMitoWorkReport;

/** Exact conservation, explicit baseline binding and missing-versus-zero evidence.
 * @author Brian, Ganyu */
public final class NcrnaMitoWorkReportTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir(),"Fixture directory must be fresh");check(new File(root+"/baseline").mkdir() && new File(root+"/C").mkdir(),"Separate source roots required");
		final String base="ncRNA Alignments (16S):\t 10\n",candidate=log("13","10","3","7");
		write(root,"baseline/stderr",base);write(root,"C/stderr",candidate);run(root,bind(root,sha,"ORIGINAL_ONLY"),"good");
		String result=read(root+"/good.tiers.tsv");check(result.contains("baseline\tP1\tORIGINAL_ONLY\t10\t10\t0\t0\tBASELINE_DERIVED") && result.contains("C\tP1\tCOVERAGE\t13\t10\t3\t7\tMEASURED"),"Only the explicit baseline derives original work and zero mitochondrial work");
		write(root,"C/stderr",base);run(root,bind(root,sha,"ORIGINAL_ONLY"),"missing");check(read(root+"/missing.tiers.tsv").contains("C\tP1\tCOVERAGE\t10\tNA\tNA\tNA\tUNAVAILABLE"),"Candidate missing counters cannot become zero");
		write(root,"baseline/stderr","No counters emitted\n");run(root,bind(root,sha,"ORIGINAL_ONLY"),"missing_base");check(read(root+"/missing_base.tiers.tsv").contains("baseline\tP1\tORIGINAL_ONLY\tNA\tNA\tNA\tNA\tUNAVAILABLE"),"Baseline derivation still requires a measured family total");
		write(root,"baseline/stderr",base);write(root,"C/stderr",log("0","0","0","0"));run(root,bind(root,sha,"ORIGINAL_ONLY"),"zero");check(read(root+"/zero.tiers.tsv").contains("C\tP1\tCOVERAGE\t0\t0\t0\t0\tMEASURED"),"Measured zero must remain distinguishable from missing");
		write(root,"C/stderr",log("14","10","3","7"));bad(root,bind(root,sha,"ORIGINAL_ONLY"),"conservation");
		write(root,"C/stderr",base+"ncRNA Original16S Alignments:\t10\n");bad(root,bind(root,sha,"ORIGINAL_ONLY"),"partial");
		write(root,"C/stderr",candidate+"ncRNA Mito Alignments:\t3\n");bad(root,bind(root,sha,"ORIGINAL_ONLY"),"duplicate");
		write(root,"C/stderr",log("13","10","bad","7"));bad(root,bind(root,sha,"ORIGINAL_ONLY"),"nonnumeric");
		write(root,"C/stderr",log("13","10","3","-1"));bad(root,bind(root,sha,"ORIGINAL_ONLY"),"negative");
		write(root,"C/stderr",log("9223372036854775808","10","3","7"));bad(root,bind(root,sha,"ORIGINAL_ONLY"),"overflow");
		write(root,"C/stderr",candidate);bad(root,bind(root,sha,"COVERAGE"),"false_baseline_policy");
		write(root,"baseline/stderr",candidate);bad(root,bind(root,sha,"ORIGINAL_ONLY"),"baseline_mito_work");
		System.out.println("MITO_WORK_REPORT_TEST_PASS measured conservation baseline_binding missing_not_zero missing_baseline partial duplicate malformed overflow");
	}
	private static String bind(String root,String sha,String policy)throws Exception{
		final String logs=NcrnaMitoCrossLibraryTest.bind(root,sha,"baseline/stderr","C/stderr");final String[] rows=read(logs).split("\n");
		write(root,"baseline/MANIFEST.sha80",rows[0].split("\t")[0]+"\tstderr\n");write(root,"C/MANIFEST.sha80",rows[1].split("\t")[0]+"\tstderr\n");
		final String manifests=NcrnaMitoCrossLibraryTest.bind(root,sha,"baseline/MANIFEST.sha80","C/MANIFEST.sha80");final String[] seals=read(manifests).split("\n");
		write(root,"runs.tsv",NcrnaMitoWorkReport.ROSTER_HEADER+"\nbaseline\tP1\t"+root+"/baseline\t"+seals[0].split("\t")[0]+"\t"+policy+"\nC\tP1\t"+root+"/C\t"+seals[1].split("\t")[0]+"\tCOVERAGE\n");
		return NcrnaMitoCrossLibraryTest.bind(root,sha,"runs.tsv","baseline/MANIFEST.sha80","C/MANIFEST.sha80","baseline/stderr","C/stderr");
	}
	private static String log(String total,String original,String mito,String skips){return "ncRNA Alignments (16S):\t"+total+"\nncRNA Original16S Alignments:\t"+original+"\nncRNA Mito Alignments:\t"+mito+"\nncRNA Mito Prealignment Length Skips:\t"+skips+"\n";}
	private static void run(String root,String pins,String name){NcrnaMitoWorkReport.main(new String[]{"runs="+root+"/runs.tsv","verified="+pins,"outprefix="+root+"/"+name,"expectedruns=2"});}
	private static void bad(String root,String pins,String name){boolean failed=false;try{run(root,pins,name);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+"/"+name+".tiers.tsv").exists(),"Invalid accounting must fail before emitting results: "+name);}
	private static void write(String root,String file,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+file,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
