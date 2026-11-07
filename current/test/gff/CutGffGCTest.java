package test.gff;

import java.io.File;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.CutGff;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GeneCaller;
import prok.NcrnaLearningContext;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Literal whole-contig GC and unchanged extraction checks for gccontig=caller.
 * CLI and batched sequence I/O follow template.A_SampleStreamer. This finite
 * fixture retains its small input records for independent cross-file comparisons.
 * @author Brian Bushnell, Raiden
 */
public final class CutGffGCTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CutGffGCTest.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires complete saved inputs and a fresh directory for the tiny-GC control. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		String ref=null, baseline=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("ref")){ref=b;}
			else if(a.equals("baseline")){baseline=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture option: "+arg);}
		}
		check(ref!=null && baseline!=null && parser.in1!=null && parser.out1!=null,
			"Specify in=caller.fa ref=genome.fa baseline=legacy.fa out=fresh_directory");
		check(parser.in2==null && parser.maxReads<0 && parser.samplerate==1 && !parser.append,
			"The fixture compares complete unpaired populations");
		final File root=new File(parser.out1);
		check(root.mkdir(), "Tiny-GC fixture output must be fresh");
		final ObjectMap<String, Read> sources=read(ref);
		final ObjectMap<String, Read> before=read(baseline);
		final ObjectMap<String, Read> after=read(parser.in1);
		check(sources.size()==3 && before.size()==5 && after.size()==5,
			"The literal population contains three source contigs and five extracted features");
		final String[] names={"ID=plus", "ID=minus", "ID=undefined", "ID=allgc_plus", "ID=allgc_minus"};
		for(final String name:names){
			final Read a=before.get(name), b=after.get(name);
			check(a!=null && b!=null && Arrays.equals(a.bases, b.bases),
				"GC annotation must preserve each extracted sequence and its orientation: "+name);
			final int x=a.id.indexOf(" contig_gc="), y=b.id.indexOf(" contig_gc=");
			check(x>0 && y>0 && a.id.substring(0, x).equals(b.id.substring(0, y)),
				"Source identity and actual clipped flanks must be unchanged: "+name);
			final String contig=RrnaResourceIO.field(b.id, "contig");
			final Read source=sources.get(contig);
			check(source!=null, "Each feature must join its complete source contig");
			final float expected=contig.equals("fraction") ? 3f/7 : contig.equals("allgc") ? 1f : 0f;
			check(bits(GeneCaller.calcContigGC(source.bases))==bits(expected),
				"Complete-contig GC must match the literal defined-base count: "+contig);
			checkGc(b, source.bases);
		}
		check(!RrnaResourceIO.field(before.get("ID=plus").id, "contig_gc").equals(
			RrnaResourceIO.field(after.get("ID=plus").id, "contig_gc")),
			"The fractional fixture must expose the legacy four-decimal rounding");
		tinyGc(root);
		System.out.println("CUTGFF_GC_TEST_PASS whole_contig float_bits native_consumer both_strands clipped_flanks undefined tiny_gc unchanged_sequences");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Verification        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercises the real teaching parser, not a replacement header interpretation. */
	private static void checkGc(final Read extracted, final byte[] source) throws Exception{
		final String text=RrnaResourceIO.field(extracted.id, "contig_gc");
		check(text.matches("[0-9]+(\\.[0-9]+)?"),
			"GC must use decimal digits with no exponent; exact zero and one need no fractional part");
		final int start=extracted.id.indexOf("contig_gc=")+10;
		final Method parse=NcrnaLearningContext.class.getDeclaredMethod("fraction", String.class, int.class, int.class);
		parse.setAccessible(true);
		final NcrnaLearningContext context=new NcrnaLearningContext();
		context.read=extracted;
		final float actual=(Float)parse.invoke(context, extracted.id, start, start+text.length());
		check(bits(actual)==bits(GeneCaller.calcContigGC(source)),
			"Serialized training GC must equal inference on the complete source contig: "+extracted.id);
	}

	/** One GC base in 10001 bases needs precision below the legacy four-decimal contract. */
	private static void tinyGc(final File root) throws Exception{
		final byte[] bases=new byte[10001];
		Arrays.fill(bases, (byte)'A');
		bases[0]='G';
		final File fasta=new File(root, "tiny.fa"), gff=new File(root, "tiny.gff"), out=new File(root, "tiny.extracted.fa");
		try(Output output=new Output(fasta.getPath())){
			output.writer.print(new ByteBuilder(">tiny\n").append(bases).nl());
		}
		try(Output output=new Output(gff.getPath())){
			output.writer.print("tiny\tfixture\tncRNA\t2\t17\t.\t+\t.\tID=tiny source=fixture contig=tiny\n");
		}
		CutGff.main(new String[]{"in="+fasta, "gff="+gff, "out="+out, "types=ncRNA",
			"flank=2", "gccontig=caller", "t=1", "ow=f"});
		final ObjectMap<String, Read> result=read(out.getPath());
		check(result.size()==1 && result.get("ID=tiny")!=null, "The tiny-GC control must retain its only feature");
		checkGc(result.get("ID=tiny"), bases);
		check(bits(GeneCaller.calcContigGC(bases))==bits((float)(1.0/10001)),
			"The tiny-GC oracle uses one GC base among 10001 defined source bases");
	}

	/** Retains only the finite fixture population through the native sequence reader. */
	private static ObjectMap<String, Read> read(final String path){
		Shared.TRIM_READ_DESCRIPTION=false;
		final ObjectMap<String, Read> result=new ObjectMap<String, Read>(String.class, Read.class);
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA,
			null, true, true), null, true, -1, false, true, 1);
		try{
			input.start();
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(final Read r:batch){
					check(r!=null && r.mate==null && result.put(RrnaResourceIO.first(r.id), r)==null,
						"Each finite fixture record must be nonnull, unpaired and uniquely named");
				}
			}
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		check(!ReadWrite.closeStream(input), "Native fixture input must close without error: "+path);
		return result;
	}

	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static void check(final boolean condition, final String reason){
		if(!condition){throw new AssertionError(reason);}
	}
}
