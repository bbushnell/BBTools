package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import prok.NcrnaMitoCoverageReport;
import prok.NcrnaMitoLengthProfile;

/** Hand-counted overlapping coverage, copy aliases, partials and bound failures.
 * @author Brian, Ganyu */
public final class NcrnaMitoLengthProfileTest {

	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be fresh");
		for(String key:NcrnaMitoLengthProfile.BINDINGS){write(root,key+".txt",key+" exact original input\n");}
		final String[] bindingFiles=new String[NcrnaMitoLengthProfile.BINDINGS.length];for(int i=0;i<bindingFiles.length;i++){bindingFiles[i]=NcrnaMitoLengthProfile.BINDINGS[i]+".txt";}
		final String bindingPins=NcrnaMitoCrossLibraryTest.bind(root,sha,bindingFiles);final StringBuilder header=new StringBuilder("#MITO_COVERAGE_MATRIX_V1\tSHORTER_QUERY_QUANTUM\tfull\n");
		final String[] pinRows=read(bindingPins).split("\n");for(int i=0;i<pinRows.length;i++){header.append("#input\t").append(NcrnaMitoLengthProfile.BINDINGS[i]).append('\t').append(pinRows[i].split("\t")[0]).append('\n');}
		header.append(NcrnaMitoLengthProfile.MATRIX_HEADER).append('\n');
		final StringBuilder matrix=new StringBuilder(header),unique=new StringBuilder();final int[] lengths={101,203,307,509};final float[][] identities={{.8f,.7f},{.66f,.9f},{.2f,.66f},{.4f,.4f}};
		for(int i=0;i<lengths.length;i++){unique.append(">U").append(i).append(" construction=true panel=").append(i==2).append(" measured=true\n").append(NcrnaMitoCrossLibraryTest.sequence(lengths[i],17+i)).append('\n');
			for(int m=0;m<3;m++){final int n=m==1?300:200;final float id=identities[i][m==1?1:0];matrix.append("U").append(i).append('\t').append(m==0?"ssu:A":m==1?"nohit:B":"ssu:Alias").append('\t').append(lengths[i]).append('\t').append(n).append('\t').append(lengths[i]<=n?"MEMBER":"MODEL").append('\t').append(id).append('\t').append(id).append('\t').append(Math.min(lengths[i],n)).append('\n');}}
		write(root,"matrix.tsv.gz",matrix.toString());write(root,"unique.fa.gz",unique.toString());
		write(root,"models.tsv",NcrnaMitoCoverageReport.MODEL_HEADER+"\nssu:A\tmitochondria_coverage_ssu_A\t200\tssu:A\tssu\tA n=3\nnohit:B\tmitochondria_coverage_nohit_B\t300\tnohit:B\tnohit\tB n=1\nssu:Alias\tmitochondria_coverage_ssu_Alias\t200\tssu:A\tssu\tAlias n=1\n");
		write(root,"selected.fa",">mitochondria_coverage_ssu_A source=ssu:A\n"+NcrnaMitoCrossLibraryTest.sequence(200,31)+"\n>mitochondria_coverage_nohit_B source=nohit:B\n"+NcrnaMitoCrossLibraryTest.sequence(300,32)+"\n");
		final String selection=NcrnaMitoCoverageReport.STEP_HEADER+"\n1\tssu:A\t2\t3\t2\t2\n2\tnohit:B\t1\t1\t3\t1\n";write(root,"selection.tsv",selection);
		final String members=NcrnaMitoCoverageReport.MEMBER_HEADER+"\n"+member("X0","U0","ssu:A","ssu:A",.8f,false)+member("X1","U0","ssu:A","ssu:A",.8f,false)+member("X2","U1","ssu:A","nohit:B",.9f,false)+member("P0","U2","nohit:B","nohit:B",.66f,true)+member("X3","U3","NA","ssu:A",.4f,false);write(root,"members.tsv.gz",members);
		write(root,"summary.tsv",NcrnaMitoCoverageReport.SUMMARY_HEADER+"\n5\t4\t3\t2\t0.66\t3\t1\t4\t1\t1\t1\tSHORTER_QUERY_COVERAGE_NOT_FULL_MODEL_CALLER_VERIFICATION\tselect\t2\t0\t12\n");
		final ArrayList<String> files=new ArrayList<String>(Arrays.asList(bindingFiles));for(String name:INPUT_NAMES){files.add(name);}
		final String validPins=pins(root,sha,files);run(root,validPins,"good");
		check(read(root+"/good.bounds.tsv").equals("model\tmodel_length\tmin_length\tmax_length\nmitochondria_coverage_ssu_A\t200\t80\t244\nmitochondria_coverage_nohit_B\t300\t80\t369\n"),"Bounds must bind model lengths in selected FASTA order and use exact floor(.8*min) and ceil(1.2*max) over ALL covered members");
		final String profile=read(root+"/good.profile.tsv");check(profile.contains("mitochondria_coverage_ssu_A\tssu:A\t200\t3\t2\t101\t203\t80\t244\t0.66\tALL_COVERED_MEMBERS"),"Model A retains two unique sequences and three copies, including the inclusive .66 hit");
		check(profile.contains("mitochondria_coverage_nohit_B\tnohit:B\t300\t4\t3\t101\t307\t80\t369\t0.66\tALL_COVERED_MEMBERS"),"B must include A-assigned members and the partial panel member; newly assigned only would give min307 and lose three copies");
		check(read(root+"/good.summary.tsv").contains("2\t5\t4\t12\t0.66\t0\tALL_COVERED_MEMBERS"),"Unique/copy universe and zero new alignments remain explicit");
		check(!new File(root+"/good.window_summary.tsv").exists() && !new File(root+"/good.window_lineages.tsv").exists(),"Omitted maxwindow must not create optional reports");
		run(root,validPins,"window180",180);
		for(String suffix:new String[]{".bounds.tsv",".profile.tsv",".summary.tsv",".inputs.tsv"}){check(read(root+"/good"+suffix).equals(read(root+"/window180"+suffix)),"Window diagnostics must leave accepted length resources and existing reports unchanged");}
		final String window=read(root+"/window180.window_summary.tsv");
		check(window.contains("180\t0.66\tINELIGIBLE_MODEL\tmitochondria_coverage_nohit_B\t300\t198\t2\t1\t4\t3\t5\t4"),"B alone exceeds the maximum window under its native .66 model-length floor");
		check(window.contains("180\t0.66\tINELIGIBLE_UNION\tALL\tNA\tNA\t2\t1\t4\t3") && window.contains("180\t0.66\tINELIGIBLE_ONLY\tALL\tNA\tNA\t2\t1\t1\t1"),"Ineligible-model coverage overlaps feasible coverage; only U2 lacks a feasible selected-model fit");
		final String lineages=read(root+"/window180.window_lineages.tsv");
		check(lineages.contains("180\t0.66\tINELIGIBLE_UNION\tALL\tA\t3\t3\tfalse") && lineages.contains("180\t0.66\tINELIGIBLE_UNION\tALL\tB\t1\t1\tfalse"),"A shared unique across lineages has two incidences but remains one unique in the global union");
		check(lineages.contains("180\t0.66\tINELIGIBLE_ONLY\tALL\tA\t1\t1\tfalse"),"Partial panel members remain in the no-feasible-coverage subset");
		run(root,validPins,"window131",131);final String both=read(root+"/window131.window_summary.tsv");
		check(both.contains("131\t0.66\tINELIGIBLE_UNION\tALL\tNA\tNA\t2\t2\t4\t3") && both.contains("131\t0.66\tINELIGIBLE_ONLY\tALL\tNA\tNA\t2\t2\t4\t3"),"Two ineligible models must union overlapping members instead of summing seven covered copies");
		run(root,validPins,"window198",198);check(read(root+"/window198.window_summary.tsv").contains("198\t0.66\tINELIGIBLE_UNION\tALL\tNA\tNA\t2\t0\t0\t0"),"A model whose native support floor equals maxwindow is feasible");
		check(NcrnaMitoLengthProfile.lowerBound(105)==84 && NcrnaMitoLengthProfile.upperBound(105)==126,"Integer multiples must not acquire floating-point off-by-one bounds");
		boolean overflow=false;try{NcrnaMitoLengthProfile.upperBound(Integer.MAX_VALUE);}catch(IllegalArgumentException e){overflow=true;}check(overflow,"Bound overflow must fail instead of wrapping");
		write(root,"matrix.tsv.gz",matrix.toString().replace("U0\tssu:A\t101", "U0\tssu:A\t102"));bad(root,pins(root,sha,files),"wrong_length");
		write(root,"matrix.tsv.gz",matrix.substring(0,matrix.lastIndexOf("U3\tssu:Alias")));bad(root,pins(root,sha,files),"missing_pair");
		write(root,"matrix.tsv.gz",matrix.toString().replace("U0\tssu:Alias\t101\t200\tMEMBER\t0.8", "U0\tssu:Alias\t101\t200\tMEMBER\t0.7"));bad(root,pins(root,sha,files),"alias_disagrees");
		write(root,"matrix.tsv.gz",matrix.toString().replace("#input\tcorpus\t","#input\twrong_key\t"));bad(root,pins(root,sha,files),"wrong_binding");
		write(root,"matrix.tsv.gz",matrix.toString());write(root,"members.tsv.gz",members+member("X0","U0","ssu:A","ssu:A",.8f,false));bad(root,pins(root,sha,files),"duplicate_copy");
		write(root,"members.tsv.gz",members.replace("X2\tU1", "X2\tU9"));bad(root,pins(root,sha,files),"missing_unique");
		write(root,"members.tsv.gz",members);write(root,"selection.tsv",selection.replace("2\tnohit:B", "2\tssu:Alias"));bad(root,pins(root,sha,files),"selected_alias");
		write(root,"selection.tsv",selection);write(root,"matrix.tsv.gz",matrix.toString().replace("U2\tnohit:B\t307\t300\tMODEL\t0.66", "U2\tnohit:B\t307\t300\tMODEL\t0.65"));bad(root,pins(root,sha,files),"assignment_disagrees");
		System.out.println("MITO_LENGTH_PROFILE_TEST_PASS all_covered_not_assigned duplicate_copies partial_panel inclusive_cut exact_rounding overflow bindings matrix_conservation alias_consistency window_overlap lineage_incidence feasible_alternatives optional_output_parity");
	}
	private static String member(String id,String uid,String assigned,String best,float score,boolean partial){final String[] p=new String[19];Arrays.fill(p,"NA");p[0]=partial?"PANEL_INCLUDED":"CONSTRUCTION";p[1]=id;p[2]=uid;p[3]="4";p[4]=id.equals("X1")?"B":"A";p[5]=""+partial;p[6]="false";p[7]=partial?"PARTIAL_INTERNAL":"ssu";p[9]=assigned;p[10]=best;p[11]=""+score;p[12]=""+score;p[13]="1";p[14]="1";p[15]=best;p[16]=""+score;p[17]=assigned.equals("NA")?"UNCOVERED_AT_CUT":"COVERED_SHORTER_QUERY";p[18]="fixture";return String.join("\t",p)+"\n";}
	private static void run(String root,String pins,String name){run(root,pins,name,-1);}
	private static void run(String root,String pins,String name,int window){final ArrayList<String> args=new ArrayList<String>();for(String key:NcrnaMitoLengthProfile.BINDINGS){args.add(key+"="+root+"/"+key+".txt");}for(int i=0;i<INPUT_KEYS.length;i++){args.add(INPUT_KEYS[i]+"="+root+"/"+INPUT_NAMES[i]);}args.add("verified="+pins);args.add("outprefix="+root+"/"+name);args.add("cut=.66");args.add("expectedselected=2");if(window>0){args.add("maxwindow="+window);}NcrnaMitoLengthProfile.main(args.toArray(new String[0]));}
	private static void bad(String root,String pins,String name){boolean failed=false;try{run(root,pins,name);}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+"/"+name+".bounds.tsv").exists(),"Bad evidence must fail before resource emission: "+name);}
	private static String pins(String root,String sha,ArrayList<String> files)throws Exception{return NcrnaMitoCrossLibraryTest.bind(root,sha,files.toArray(new String[0]));}
	private static void write(String root,String name,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+name,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] INPUT_KEYS={"matrix","unique","models","members","selection","summary","selected"};
	private static final String[] INPUT_NAMES={"matrix.tsv.gz","unique.fa.gz","models.tsv","members.tsv.gz","selection.tsv","summary.tsv","selected.fa"};
}
