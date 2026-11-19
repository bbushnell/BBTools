package test.prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import ddl.FamilySketchBuilder;
import dna.AminoAcid;
import fileIO.ByteFile;
import prot.Blosum62;
import prot.HbmPositionModel;
import prot.ProteinSequence;
import prot.ProteinSubstitutionMatrix;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.DdlClust;
import prot.cluster.Family;
import prot.cluster.TailAlignmentDecision;
import prot.cluster.TailAlignmentGate;

/** Native CLI controls against retained parent outputs and a literal changed-trace
 * witness. No biological population or calibration is consumed.
 * @author Brian, Qiqi
 */
public final class ClusteringMatrixTest{

	public static void main(String[] args)throws Exception{
		check(args.length==2 || args.length==4,"Supply fresh directory and parent, or candidate plus parent directory and empirical matrix");
		final boolean parent=args[1].equals("parent");
		check(parent ? args.length==2 : args[1].equals("candidate") && args.length==4,"Explicit fixture mode required");
		final Path root=Paths.get(args[0]);
		check(!Files.exists(root),"Never overwrite retained native fixture evidence");Files.createDirectories(root);
		final String a=repeat('A',100),w=repeat('W',100),c=repeat('C',90);
		final Path input=root.resolve("parity.faa"),seeds=root.resolve("parity.pddl");
		write(input,">a0\n"+a+"\n>a1\n"+a+"\n>a2\n"+a+"\n>w\n"+w+"\n>short\n"+c+"\n");
		seedCodec(seeds,new String[]{a,w,repeat('C',100)});
		// Retain the historical artifact names, but bind the legacy control explicitly.
		greedy(root.resolve("default_greedy"),input,"blosum62",null,100,70,1,false,null);
		cake(root.resolve("default_cake"),input,seeds,"blosum62",null,1);
		if(parent){
			finish(root,"CLUSTERING_MATRIX_PARENT_PASS");return;
		}
		final Path old=Paths.get(args[2]);
		check(Files.readAllLines(old.resolve("PASS"),StandardCharsets.US_ASCII).contains("CLUSTERING_MATRIX_PARENT_PASS"),"Parent evidence must precede candidate comparisons");
		parity(old.resolve("default_greedy"),root.resolve("default_greedy"),false);
		parity(old.resolve("default_cake"),root.resolve("default_cake"),true);
		final ProteinSubstitutionMatrix scaled=ProteinSubstitutionMatrix.fromScores("custom40",4,16,scores());
		final int[][] changedScores=scores();
		final int ai=AminoAcid.acidToNumberExtended['A'],ri=AminoAcid.acidToNumberExtended['R'];
		changedScores[ai][ri]=changedScores[ri][ai]=16;
		final ProteinSubstitutionMatrix changed=ProteinSubstitutionMatrix.fromScores("custom60",4,16,changedScores);
		final Path scaledFile=root.resolve("scaled.tsv"),changedFile=root.resolve("changed.tsv");
		write(scaledFile,resource(scaled));write(changedFile,resource(changed));
		for(String name:new String[]{"blosum62","custom40"}){
			final Path matrix=name.equals("custom40") ? scaledFile : null;
			greedy(root.resolve(name+"_greedy"),input,name,matrix,100,70,1,false,null);
			parity(old.resolve("default_greedy"),root.resolve(name+"_greedy"),false);
			for(int threads:new int[]{1,4}){
				final Path out=root.resolve(name+"_cake_t"+threads);
				cake(out,input,seeds,name,matrix,threads);
				parity(old.resolve("default_cake"),out,true);
				if(matrix!=null){binding(out.resolve("final/stats.tsv"),scaled);}
			}
			if(matrix!=null){binding(root.resolve(name+"_greedy/stats.tsv"),scaled);}
		}
		traceOracle(root,changed);
		final Path witness=root.resolve("witness.faa"),oneSeed=root.resolve("witness.pddl");
		final StringBuilder fasta=new StringBuilder();
		for(int i=0;i<12;i++){fasta.append('>').append('q').append(i).append('\n').append(i<6 ? "AAAA" : "RAAA").append('\n');}
		write(witness,fasta.toString());seedCodec(oneSeed,new String[]{"AAAA"});
		for(boolean concurrent:new boolean[]{false,true}){
			for(String name:new String[]{"blosum62","custom60"}){
				final Path out=root.resolve("witness_"+name+"_"+concurrent);
				greedy(out,witness,name,name.equals("custom60") ? changedFile : null,4,4,.75,concurrent,null);
				final ClusterSketchCodec.Loaded models=ClusterSketchCodec.read(out.resolve("families.pddl.gz").toString());
				check(models.families.length==(name.equals("custom60") ? 1 : 2),"Trace witness requires custom one family versus BLOSUM two; mode="+concurrent);
				check(assigned(out.resolve("members.tsv.gz"))==12,"Greedy must retain all twelve witness inputs");
			}
		}
		for(String name:new String[]{"blosum62","custom60"}){
			final Path out=root.resolve("frozen_"+name);
			greedy(out,witness,name,name.equals("custom60") ? changedFile : null,4,4,.75,false,oneSeed);
			check(assigned(out.resolve("members.tsv.gz"))==(name.equals("custom60") ? 12 : 6),"Frozen AAAA seed must accept RAAA only under changed optimum");
		}
		final Path empiricalFile=Paths.get(args[3]).toAbsolutePath().normalize();
		check(empiricalFile.toString().equals(ProteinSubstitutionMatrix.defaultResourcePath()),
			"The default resolver must select the exact installed POSSUM55 resource");
		retain(empiricalFile);
		final ProteinSubstitutionMatrix empirical=ProteinSubstitutionMatrix.load("POSSUM55",empiricalFile.toString());
		check("POSSUM55".equals(ProteinSubstitutionMatrix.DEFAULT_NAME),"New alignment work defaults to POSSUM55");
		final ProteinSubstitutionMatrix implicit=ProteinSubstitutionMatrix.defaultMatrix();
		final ProteinSubstitutionMatrix namedDefault=ProteinSubstitutionMatrix.load("POSSUM55",null);
		check(empirical.sha80.equals("21a59255f7b0dfa77945") && empirical.scale==4 && empirical.linearGap==16,
			"The empirical control is the accepted POSSUM55 binding, scale and gap");
		// Literal resource witnesses distinguish POSSUM55 from the historical custom40 scores.
		check(empirical.score((byte)ai,(byte)ai)==17 && empirical.score((byte)ai,(byte)ri)==-7 &&
			empirical.score((byte)ri,(byte)ai)==-7 && empirical.score((byte)ai,(byte)AminoAcid.acidToNumberExtended['N'])==-13,
			"Installed POSSUM55 must retain A/A=17, A/R=R/A=-7 and A/N=-13, not merely a renamed default");
		for(ProteinSubstitutionMatrix matrix:new ProteinSubstitutionMatrix[]{implicit,namedDefault}){
			check(matrix.name.equals(empirical.name) && matrix.sha80.equals(empirical.sha80) &&
				matrix.scale==empirical.scale && matrix.linearGap==empirical.linearGap,"Central implicit and named defaults bind the accepted matrix");
		}
		greedy(root.resolve("empirical_greedy"),input,"POSSUM55",empiricalFile,100,70,1,false,null);
		cake(root.resolve("empirical_cake"),input,seeds,"POSSUM55",empiricalFile,2);
		binding(root.resolve("empirical_greedy/stats.tsv"),empirical);
		binding(root.resolve("empirical_cake/final/stats.tsv"),empirical);
		for(boolean named:new boolean[]{false,true}){
			final String label=named ? "named_default" : "implicit_default";
			greedy(root.resolve(label+"_greedy"),input,named ? "POSSUM55" : null,null,100,70,1,false,null);
			cake(root.resolve(label+"_cake"),input,seeds,named ? "POSSUM55" : null,null,2);
			binding(root.resolve(label+"_greedy/stats.tsv"),empirical);
			binding(root.resolve(label+"_cake/final/stats.tsv"),empirical);
			parity(root.resolve("empirical_greedy"),root.resolve(label+"_greedy"),false);
			parity(root.resolve("empirical_cake"),root.resolve(label+"_cake"),true);
			for(String recipe:new String[]{"run.config","layer1/run.config","layer2/run.config"}){
				check(setting(root.resolve(label+"_cake").resolve(recipe),"matrixfile","=").equals(empiricalFile.toString()),
					"Default Cake replay recipes must retain the exact resolved matrix input");
			}
		}
		implicitSelectorAliases(root,input,empiricalFile);
		rejections(root,input,seeds,scaledFile);
		finish(root,"CLUSTERING_MATRIX_TEST_PASS parent_parity=true scaled_flow=true changed_trace=true serial=true concurrent=true frozen=true cake=true empirical=true POSSUM55_default=true explicit_blosum=true implicit_selector_file_aliases=true rejected=true");
	}

	/** Checks implicit matrix selection with an explicit fresh input copy. Installed
	 * default-path resolution is checked separately; no installed resource is writable here.
	 */
	private static void implicitSelectorAliases(Path root,Path input,Path empirical)throws Exception{
		final byte[] original=Files.readAllBytes(empirical);
		final Path copy=root.resolve("alias_matrix.tsv").toAbsolutePath().normalize();
		final Path link=root.resolve("alias_matrix_hardlink.tsv").toAbsolutePath().normalize();
		Files.write(copy,original,StandardOpenOption.CREATE_NEW);Files.createLink(link,copy);
		final StringBuilder errors=new StringBuilder();
		for(int i=0;i<2;i++){
			final Path target=i==0 ? copy : link,out=root.resolve("alias_rejected_"+i);
			final ArrayList<String> args=greedyArgs(out,input,100,70,1);
			args.add("matrixfile="+copy);
			for(int j=0;j<args.size();j++){
				if(args.get(j).startsWith("stats=")){args.set(j,"stats="+target);}
				if(args.get(j).equals("ow=f")){args.set(j,"ow=t");}
			}
			write(root.resolve("alias_rejected_"+i+".argv"),String.join("\n",args)+"\n");
			boolean rejected=false;
			try{DdlClust.main(args.toArray(new String[args.size()]));}
			catch(IllegalArgumentException failure){
				final String diagnostic=failure.getMessage();
				check(diagnostic!=null && diagnostic.toLowerCase(Locale.ROOT).contains("alias") && diagnostic.contains(target.toString()),
					"Only the exact matrix-file alias diagnostic satisfies this rejection: "+failure);
				errors.append(i).append('\t').append(diagnostic.replace('\n',' ')).append('\n');rejected=true;
			}
			check(rejected && !Files.exists(out),"Implicit-selector matrix-file aliases must fail before creating outputs");
			check(Arrays.equals(original,Files.readAllBytes(copy)) && Arrays.equals(original,Files.readAllBytes(link)),
				"Neither exact-path nor hard-link output may change the fixture matrix");
		}
		check(Arrays.equals(original,Files.readAllBytes(empirical)),"Installed POSSUM55 resource remains unchanged");
		retain(copy);retain(link);write(root.resolve("implicit_selector_file_alias_rejections.tsv"),errors.toString());
	}

	/** Published BLOSUM A/R=-1 clips the first pair; changed +16 consumes it. */
	private static void traceOracle(Path root,ProteinSubstitutionMatrix changed)throws Exception{
		final ProteinSequence a=new ProteinSequence("a","AAAA"),r=new ProteinSequence("r","RAAA");
		final TailAlignmentGate old=gate(ProteinSubstitutionMatrix.BLOSUM62),custom=gate(changed);
		final StringBuilder receipt=new StringBuilder();
		for(boolean reverse:new boolean[]{false,true}){
			final byte[] query=reverse ? a.enc : r.enc,reference=reverse ? r.enc : a.enc;
			final TailAlignmentDecision x=old.decide(query,reference),y=custom.decide(query,reference);
			check(x.reason==TailAlignmentDecision.Reason.QUERY_SPAN_SHORT && x.queryStart==1 && x.querySpan==3 && x.score==12*HbmPositionModel.SCALE,
				"BLOSUM literal optimum is three A pairs, excluded by the span-four gate");
			check(y.accepted() && y.queryStart==0 && y.referenceStart==0 && y.querySpan==4 && y.referenceSpan==4 &&
				y.identities==3 && y.paired==4 && y.identity==.75 && y.score==64*HbmPositionModel.SCALE,
				"Changed A/R score must produce four paired columns before any clustering count oracle is used");
			receipt.append(reverse).append('\t').append(x).append('\t').append(y).append('\n');
		}
		write(root.resolve("trace_oracle.tsv"),receipt.toString());
	}
	private static TailAlignmentGate gate(ProteinSubstitutionMatrix matrix){
		return new TailAlignmentGate(TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT,TailAlignmentGate.Denominator.ALIGNED_QUERY,100,4,.75,0,null,matrix);
	}
	private static ArrayList<String> greedyArgs(Path out,Path input,int minLength,int span,double identity){
		check(out!=null && input!=null,"Explicit fixture paths required");
		return new ArrayList<String>(Arrays.asList("in="+input,"memberranking=f","t=1","seed=1701","k=2","width=64","sketchseed=7",
			"positional=f","partitions=0","overlap=0","blacklist=null","minlen="+minLength,"minspan="+span,"maxalignlen=100",
			"minid="+identity,"minblosumratio=0","pinright=t","trimleft=t","consensustrimleft=t","denominator=aligned",
			"mincluster=1","passes=1","consensuspasses=2","merge=f","indexed=t","greedyhits=0","passhits=0",
			"maxhydrodelta=9","greedymaxhydrodelta=9","out=null","members="+out.resolve("members.tsv.gz"),
			"ddlout="+out.resolve("families.pddl.gz"),"consensus="+out.resolve("consensus.faa.gz"),"stats="+out.resolve("stats.tsv"),
			"quality="+out.resolve("quality"),"refinement=null","checkpoints=null","hist=null","pairhist=null","ow=f"));
	}
	private static void greedy(Path out,Path input,String name,Path matrix,int minLength,int span,double identity,boolean concurrent,Path seeds)throws Exception{
		Files.createDirectory(out);
		final ArrayList<String> args=greedyArgs(out,input,minLength,span,identity);select(args,name,matrix);
		if(concurrent){args.addAll(Arrays.asList("concurrentgreedy=t","greedythreads=1","snapshotfraction=.01"));}
		if(seeds!=null){args.add("seeds="+seeds);}
		invoke(out,args);
		for(String leaf:new String[]{"members.tsv.gz","families.pddl.gz","consensus.faa.gz","stats.tsv","quality.cohesion.tsv"}){retain(out.resolve(leaf));}
	}
	private static ArrayList<String> cakeArgs(Path out,Path input,Path seeds,int threads){
		check(threads>0,"Native test threads must be positive");
		return new ArrayList<String>(Arrays.asList("mode=cake","memberranking=f","in="+input,"codecs="+seeds,"outdir="+out,
			"familybudget=1","initialfraction=3","passes=2","t="+threads,"k=2","width=64","sketchseed=7",
			"positional=f","partitions=0","overlap=0","blacklist=null","minlen=100","minspan=70","maxalignlen=100",
			"minid=1","minblosumratio=0","mergeminid=1","mergeminratio=0","trimleft=f","consensustrimleft=f",
			"passhits=0","greedyhits=0","maxhydrodelta=9","greedymaxhydrodelta=9","fractions=.25,1",
			"sizelambda=.125","belowmode=constant","minchildsize=1","fissiontopfraction=1","fusionbottomfraction=1"));
	}
	private static void cake(Path out,Path input,Path seeds,String name,Path matrix,int threads)throws Exception{
		final ArrayList<String> args=cakeArgs(out,input,seeds,threads);select(args,name,matrix);
		// Cake owns output-directory creation; retain argv beside it before invoking.
		write(out.resolveSibling(out.getFileName()+".argv"),String.join("\n",args)+"\n");
		DdlClust.main(args.toArray(new String[args.size()]));
		retain(out.resolve("layers.tsv"));retain(out.resolve("run.config"));
		for(int layer=1;layer<=2;layer++){
			final Path stage=out.resolve("layer"+layer);
			final List<String> waves=Files.readAllLines(stage.resolve("waves.tsv"),StandardCharsets.US_ASCII);
			boolean fusion=false;for(String row:waves){if(row.startsWith("fusion\t")){fusion=true;}}
			check(fusion,"Both resident layers must execute the fusion phase");
			check(Integer.parseInt(setting(stage.resolve("state.tsv"),"fission_waves","\t"))>0,"Tiny resident fixture must enter fission waves");
			for(String leaf:new String[]{"waves.tsv","state.tsv","run.config"}){retain(stage.resolve(leaf));}
		}
		for(String leaf:new String[]{"members.tsv.gz","families.pddl.gz","consensus.faa.gz","stats.tsv","quality.cohesion.tsv"}){retain(out.resolve("final").resolve(leaf));}
		final ClusterSketchCodec.Loaded result=ClusterSketchCodec.read(out.resolve("final/families.pddl.gz").toString());
		check(result.families.length==1 && result.families[0].members==3,"Final top-one reassignment must retain the three A sequences");
	}
	private static void invoke(Path out,ArrayList<String> args)throws Exception{
		write(out.resolve("argv.txt"),String.join("\n",args)+"\n");DdlClust.main(args.toArray(new String[args.size()]));
	}
	private static void select(ArrayList<String> args,String name,Path matrix){
		check(matrix==null || name!=null,"Resource cannot have an implicit selector in the fixture");
		if(name!=null){args.add("matrix="+name);}if(matrix!=null){args.add("matrixfile="+matrix);}
	}
	private static void rejections(Path root,Path input,Path seeds,Path scaled)throws Exception{
		final Path malformed=root.resolve("malformed.tsv");write(malformed,"not a matrix\n");
		final String[][] variants={{"custom60",null},{"custom40",root.resolve("absent.tsv").toString()},
			{"custom40",malformed.toString()},{"custom60",scaled.toString()},{"custom40",root.toString()}};
		final StringBuilder errors=new StringBuilder();
		for(boolean cake:new boolean[]{false,true}){
			for(int i=0;i<variants.length;i++){
				final Path out=root.resolve("rejected_"+cake+"_"+i);
				final ArrayList<String> args=cake ? cakeArgs(out,input,seeds,1) : greedyArgs(out,input,100,70,1);
				select(args,variants[i][0],variants[i][1]==null ? null : Paths.get(variants[i][1]));
				boolean rejected=false;
				try{DdlClust.main(args.toArray(new String[args.size()]));}
				catch(IllegalArgumentException failure){
					check(failure.getMessage()!=null,"Invalid resource diagnostic must explain the rejection");
					final String diagnostic=failure.getMessage().toLowerCase(Locale.ROOT);
					check(diagnostic.contains("matrix") || diagnostic.contains("substitution"),"Unrelated path/argument failure cannot satisfy a matrix rejection: "+failure);
					errors.append(cake).append('\t').append(i).append('\t').append(failure.getClass().getName()).append('\t')
						.append(failure.getMessage().replace('\n',' ')).append('\n');rejected=true;
				}
				check(rejected && !Files.exists(out),"Missing/malformed/mismatched matrix must fail before any output directory exists");
			}
		}
		write(root.resolve("rejections.tsv"),errors.toString());
	}
	private static void parity(Path first,Path second,boolean cake)throws Exception{
		if(cake){same(first.resolve("layers.tsv"),second.resolve("layers.tsv"));first=first.resolve("final");second=second.resolve("final");}
		for(String leaf:new String[]{"members.tsv.gz","families.pddl.gz","consensus.faa.gz","quality.cohesion.tsv"}){same(first.resolve(leaf),second.resolve(leaf));}
	}
	private static void same(Path first,Path second){
		final ByteFile a=ByteFile.makeByteFile(first.toString(),true),b=ByteFile.makeByteFile(second.toString(),true);
		try{
			while(true){final byte[] x=a.nextLine(),y=b.nextLine();check(Arrays.equals(x,y),"Decoded result parity: "+first+" versus "+second);if(x==null){break;}}
		}finally{final boolean x=a.close(),y=b.close();check(!x && !y,"Both decoded evidence streams must close cleanly");}
		retain(first);retain(second);
	}
	private static int assigned(Path members){
		final ByteFile input=ByteFile.makeByteFile(members.toString(),true);int count=0;
		try{for(byte[] row=input.nextLine();row!=null;row=input.nextLine()){
			if(new String(row,StandardCharsets.US_ASCII).endsWith("\tassigned")){count++;}
		}}finally{check(!input.close(),"Membership stream must close cleanly");}return count;
	}
	private static void binding(Path stats,ProteinSubstitutionMatrix matrix)throws Exception{
		check(setting(stats,"alignment_matrix","\t").equals(matrix.name) && setting(stats,"alignment_matrix_sha80","\t").equals(matrix.sha80) &&
			setting(stats,"alignment_matrix_scale","\t").equals(Integer.toString(matrix.scale)) &&
			setting(stats,"alignment_matrix_gap","\t").equals(Integer.toString(matrix.linearGap)),"Actual runtime statistics must identify the selected immutable matrix and its units");
	}
	private static String setting(Path file,String key,String delimiter)throws Exception{
		String value=null;for(String line:Files.readAllLines(file,StandardCharsets.US_ASCII)){
			if(line.startsWith(key+delimiter)){check(value==null,"Duplicate evidence field: "+key);value=line.substring(key.length()+delimiter.length());}
		}check(value!=null,"Missing evidence field: "+key+" in "+file);return value;
	}
	private static void seedCodec(Path path,String[] texts){
		check(texts.length>0,"At least one literal seed is required");
		final ProteinSequence[] seqs=new ProteinSequence[texts.length];final Family[] families=new Family[texts.length];
		final FamilySketchBuilder builder=new FamilySketchBuilder(64,2,7,false,0,0,null);
		for(int i=0;i<texts.length;i++){seqs[i]=new ProteinSequence("seed"+i,texts[i]);families[i]=new Family(7+i*12,"seed"+i,i,seqs[i].enc);}
		for(Family family:families){family.refreshStats(seqs,builder);}
		ClusterSketchCodec.write(path.toString(),false,families,2,64,7,false,0,0);retain(path);
	}
	private static int[][] scores(){
		final int[][] result=new int[20][20];for(int a=0;a<20;a++){for(int b=0;b<20;b++){result[a][b]=4*Blosum62.score((byte)a,(byte)b);}}return result;
	}
	private static String resource(ProteinSubstitutionMatrix matrix){
		check(matrix!=null,"Explicit matrix required");
		final StringBuilder out=new StringBuilder("#bbgap_substitution_matrix\t1\nname\t").append(matrix.name)
			.append("\nscale\t").append(matrix.scale).append("\ngap\t").append(matrix.linearGap).append('\n');
		for(char a:"ARNDCQEGHILKMFPSTWYV".toCharArray()){for(char b:"ARNDCQEGHILKMFPSTWYV".toCharArray()){
			out.append("score\t").append(a).append('\t').append(b).append('\t').append(matrix.score((byte)AminoAcid.acidToNumberExtended[a],(byte)AminoAcid.acidToNumberExtended[b])).append('\n');
		}}return out.toString();
	}
	private static String repeat(char c,int n){check(n>0,"Positive motif length required");final char[] chars=new char[n];Arrays.fill(chars,c);return new String(chars);}
	private static void write(Path path,String text)throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII),StandardOpenOption.CREATE_NEW);retain(path);}
	private static void retain(Path path){check(Files.isRegularFile(path),"Expected retained evidence: "+path);final String name=path.toAbsolutePath().normalize().toString();if(!evidence.contains(name)){evidence.add(name);}}
	private static void finish(Path root,String marker)throws Exception{
		write(root.resolve("PASS"),marker+"\n");Files.write(root.resolve("evidence.paths"),evidence,StandardCharsets.US_ASCII,StandardOpenOption.CREATE_NEW);System.out.println(marker);
	}
	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}}
	private static final ArrayList<String> evidence=new ArrayList<String>();
	private ClusteringMatrixTest(){}
}
