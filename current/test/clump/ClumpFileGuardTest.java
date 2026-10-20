package test.clump;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import clump.KmerSort2;
import clump.KmerSort3;
import clump.KmerSplit;
import clump.KmerReduce;
import clump.Clumpify;
import fileIO.ReadWrite;
import shared.Shared;
import stream.FASTQ;

/** Constructor-only regression cases for clump file alias guards.
 * Exercises paired and percent-expanded names without processing reads.
 * @author Jean, Shinobu
 */
public class ClumpFileGuardTest{

	/** Runs baseline or repaired expectations in a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		check(args.length==2 || (args.length==3 &&
				(args[2].equals("out2") || args[2].equals("orchestrator"))),
				"Expected baseline|fixed, a fresh fixture directory, and optional out2|orchestrator");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final StaticState state=new StaticState();
		try{
			state.setDeterministic();
			final Path root=Paths.get(args[1]);
			Files.createDirectory(root);
			System.out.println("case\taccepted\tinputs_unchanged\tdistinct_outputs_absent");
			if(args.length==3){
				if(args[2].equals("orchestrator")){
					for(int mode=0; mode<10; mode++){checkOrchestrator(root, mode, fixed);}
				}else{
					for(int mode=0; mode<5; mode++){checkMateOutput(root, mode, fixed);}
				}
				return;
			}
			checkPaired(root, "kmersort2-paired-alias", false, true, fixed);
			checkPaired(root, "kmersort3-paired-alias", true, true, fixed);
			checkExpanded(root, "kmersort2-expanded-alias", false, true, fixed);
			checkExpanded(root, "kmersort3-expanded-alias", true, true, fixed);
			checkPaired(root, "kmersort2-paired-distinct", false, false, fixed);
			checkPaired(root, "kmersort3-paired-distinct", true, false, fixed);
			checkExpanded(root, "kmersort2-expanded-distinct", false, false, fixed);
			checkExpanded(root, "kmersort3-expanded-distinct", true, false, fixed);
			for(boolean splitter : new boolean[] {false, true}){
				for(int collision=-1; collision<2; collision++){
					checkUtility(root, splitter, collision, fixed);
				}
			}
		}finally{
			state.restore();
		}
	}

	/** Checks original inputs against final destinations without starting any processing. */
	private static void checkOrchestrator(Path root, int mode, boolean fixed) throws Exception{
		final String[] names={"out1-in1", "out1-in2", "out2-in1", "out2-in2",
				"same-outputs", "expanded-hash", "distinct-grouped", "distinct-single", "no-output", "null-options"};
		final Path dir=root.resolve(names[mode]);
		Files.createDirectory(dir);
		final Path left=writeInput(dir, "input1.fq"), right=writeInput(dir, "input2.fq");
		final Path one=dir.resolve("one.fq"), two=dir.resolve("two.fq");
		String in1=left.toString(), in2=right.toString(), out1=one.toString(), out2=two.toString();
		if(mode==0){out1=in1;}
		else if(mode==1){out1=in2;}
		else if(mode==2){out2=in1;}
		else if(mode==3){out2=in2;}
		else if(mode==4){out2=out1;}
		else if(mode==5){in1=out1=dir.resolve("input#.fq").toString(); in2=out2=null;}
		else if(mode==8){out1=out2=null;}
		else if(mode==9){in2=out2="null";}
		final ArrayList<String> args=new ArrayList<String>();
		args.add("in="+in1);
		if(in2!=null){args.add("in2="+in2);}
		if(out1!=null){args.add("out="+out1);}
		if(out2!=null){args.add("out2="+out2);}
		args.add("groups="+(mode==7 ? 1 : 2));
		args.add("overwrite=t");
		args.add("deleteinput=t");
		boolean accepted=true;
		try{new Clumpify(args.toArray(new String[args.size()]));}
		catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains("specified multiple times"),
					"Unexpected orchestrator constructor rejection: "+e);
			accepted=false;
		}
		check(accepted==(!fixed || mode>=6), "Unexpected orchestrator acceptance: "+names[mode]);
		checkInput(dir, "input1.fq");
		checkInput(dir, "input2.fq");
		check(!Files.exists(one) && !Files.exists(two), "Constructor created output files");
		check(!Files.exists(dir.resolve("input#.fq")), "Constructor created a literal hash filename");
		System.out.println(names[mode]+'\t'+accepted+"\ttrue\ttrue");
	}

	/** Checks the second output's required primary, overwrite and alias contracts. */
	private static void checkMateOutput(Path root, int mode, boolean fixed) throws Exception{
		final String[] names={"out2-only", "out2-existing", "out2-input-alias", "out2-distinct", "out2-append"};
		final Path dir=root.resolve(names[mode]);
		Files.createDirectory(dir);
		final Path input=writeInput(dir, "input.fq"), primary=dir.resolve("one.fq");
		final Path mate=mode==2 ? input : dir.resolve("two.fq");
		if(mode==1 || mode==4){writeInput(dir, "two.fq");}
		final String[] args={"in="+input, "out="+(mode==0 ? "null" : primary), "out2="+mate,
				"groups=1", "overwrite="+(mode==2), "append="+(mode==4)};
		boolean accepted=true;
		try{new KmerSort2(args);}
		catch(RuntimeException | AssertionError e){
			final String message=String.valueOf(e.getMessage()).toLowerCase(java.util.Locale.ROOT);
			check(mode==0 ? message.contains("out1") && message.contains("out2") :
					mode==1 ? message.contains("write") : mode==2 && message.contains("specified multiple times"),
					"Unexpected second-output rejection: "+e);
			accepted=false;
		}
		check(accepted==(!fixed || mode>=3), "Unexpected second-output acceptance: "+names[mode]);
		checkInput(dir, "input.fq");
		if(mode==1 || mode==4){checkInput(dir, "two.fq");}
		check(!Files.exists(primary), "Constructor created the primary output");
		if(mode==0 || mode==3){check(!Files.exists(mate), "Constructor created the mate output");}
		System.out.println(names[mode]+'\t'+accepted+"\ttrue\ttrue");
	}

	/** Ensures out=in2 aliasing is rejected while distinct paired files remain valid. */
	private static void checkPaired(Path root, String name, boolean kmerSort3, boolean alias,
			boolean fixed) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path in1=writeInput(dir, "left.fq"), in2=writeInput(dir, "right.fq");
		final Path out=alias ? in2 : dir.resolve("out.fq");
		final String[] args={"in1="+in1, "in2="+in2, "out="+out, "groups=1"};
		final boolean accepted=construct(args, kmerSort3);
		check(accepted!=(alias && fixed), "Unexpected paired alias result for "+name);
		checkInput(dir, "left.fq");
		checkInput(dir, "right.fq");
		final boolean outputsAbsent=alias || !Files.exists(out);
		check(outputsAbsent, "Constructor-only paired case created output: "+out);
		System.out.println(name+'\t'+accepted+'\t'+true+'\t'+outputsAbsent);
	}

	/** Ensures aliases visible only after % expansion are rejected. */
	private static void checkExpanded(Path root, String name, boolean kmerSort3, boolean alias,
			boolean fixed) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		writeInput(dir, "in_0.fq");
		writeInput(dir, "in_1.fq");
		final Path inPattern=dir.resolve("in_%.fq");
		final Path outPattern=dir.resolve(alias ? "in_%.fq" : "out_%.fq");
		final String[] args={"in="+inPattern, "out="+outPattern, "groups=2"};
		final boolean accepted=construct(args, kmerSort3);
		check(accepted!=(alias && fixed), "Unexpected expanded alias result for "+name);
		checkInput(dir, "in_0.fq");
		checkInput(dir, "in_1.fq");
		final boolean outputsAbsent=alias ||
				(!Files.exists(dir.resolve("out_0.fq")) && !Files.exists(dir.resolve("out_1.fq")));
		check(outputsAbsent, "Constructor-only expanded case created output: "+outPattern);
		System.out.println(name+'\t'+accepted+'\t'+true+'\t'+outputsAbsent);
	}

	/** Constructs one sorter without invoking process(). */
	private static boolean construct(String[] args, boolean kmerSort3){
		try{
			if(kmerSort3){new KmerSort3(0, 1, 1, args);}
			else{new KmerSort2(args);}
			return true;
		}catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains("specified multiple times"),
					"Unexpected constructor rejection: "+e);
			return false;
		}
	}

	/** Checks both input aliases and distinct destinations for the split/reduce frontends. */
	private static void checkUtility(Path root, boolean splitter, int collision, boolean fixed) throws Exception{
		final String name=(splitter ? "kmersplit" : "kmerreduce")+"-"+
				(collision<0 ? "distinct" : collision==0 ? "primary-alias" : "mate-alias");
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path in1=writeInput(dir, "left0.fq"), in2=writeInput(dir, "right0.fq");
		final String prefix=collision<0 ? "out" : collision==0 ? "left" : "right";
		final Path output=dir.resolve(prefix+(splitter ? "%" : "0")+".fq");
		final String[] args={"in1="+in1, "in2="+in2, "out="+output, "groups=2", "overwrite=t"};
		boolean accepted=true;
		try{
			if(splitter){new KmerSplit(args);}else{new KmerReduce(args);}
		}catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains("specified multiple times"),
					"Unexpected utility constructor rejection: "+e);
			accepted=false;
		}
		check(accepted!=((collision>=0) && fixed), "Unexpected utility alias result: "+name);
		checkInput(dir, "left0.fq");
		checkInput(dir, "right0.fq");
		if(collision<0){check(!Files.exists(dir.resolve("out0.fq")), "Constructor created output");}
		if(splitter){check(!Files.exists(dir.resolve(prefix+"1.fq")), "Constructor created a second group");}
		System.out.println(name+'\t'+accepted+"\ttrue\ttrue");
	}

	/** Writes one FASTQ fixture with CREATE_NEW so no existing data can be replaced. */
	private static Path writeInput(Path dir, String name) throws Exception{
		final Path path=dir.resolve(name);
		Files.write(path, FASTQ_TEXT.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		return path;
	}

	/** Verifies constructor-only cases did not overwrite a disposable input. */
	private static void checkInput(Path dir, String name) throws Exception{
		final Path path=dir.resolve(name);
		check(new String(Files.readAllBytes(path), StandardCharsets.US_ASCII).equals(FASTQ_TEXT),
				"Input changed during constructor-only alias check: "+path);
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Static state touched by sorter constructors. */
	private static class StaticState{
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final float zipThreadMult=ReadWrite.ZIP_THREAD_MULT();
		final int maxZipThreads=ReadWrite.MAX_ZIP_THREADS(), zipLevel=ReadWrite.ZIPLEVEL;
		final boolean forceInterleaved=FASTQ.FORCE_INTERLEAVED, testInterleaved=FASTQ.TEST_INTERLEAVED;

		void setDeterministic(){
			Shared.setThreads(1);
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		void restore(){
			ReadWrite.USE_PIGZ=pigz;
			ReadWrite.USE_UNPIGZ=unpigz;
			ReadWrite.setZipThreadMult(zipThreadMult);
			ReadWrite.setZipThreads(maxZipThreads);
			ReadWrite.ZIPLEVEL=zipLevel;
			FASTQ.FORCE_INTERLEAVED=forceInterleaved;
			FASTQ.TEST_INTERLEAVED=testInterleaved;
		}
	}

	private static final String FASTQ_TEXT="@alias\nACGTACGT\n+\nIIIIIIII\n";
}
