package clade;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import bin.AdjustEntropy;
import cardinality.DynamicDemiLog;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.EntropyTracker;

/** Serial byte equivalence across worker counts, Read containers and file loading.
 * @author Ganyu */
public final class CladeBuilderTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1 && args[0].startsWith("outdir="),"Fresh outdir required");
		final Path root=Files.createDirectory(Paths.get(args[0].substring(7)));
		AdjustEntropy.load(4,150);Clade.MAKE_DDLS=true;Clade.DDL_K=25;Clade.DDL_BUCKETS=32768;DynamicDemiLog.setExponent(5);
		final ArrayList<Read> reads=new ArrayList<Read>();final Random random=new Random(818);
		for(int i=0;i<48;i++){
			final byte[] bases=new byte[200+random.nextInt(8000)];
			for(int j=0;j<bases.length;j++){bases[j]=(byte)"ACGTNacgt".charAt(random.nextInt(9));}
			reads.add(new Read(bases,null,"contig"+i,i));
		}
		reads.add(new Read("NNNNNNNNNN".getBytes(StandardCharsets.US_ASCII),null,"allN",48));
		reads.add(new Read("AC".getBytes(StandardCharsets.US_ASCII),null,"short",49));
		final byte[][] before=new byte[reads.size()][];
		for(int i=0;i<reads.size();i++){before[i]=reads.get(i).bases.clone();}
		for(boolean entropy:new boolean[]{false,true}){
			CladeObject.calcCladeEntropy=entropy;
			final Clade expected=serial(reads,"fixture");
			for(int threads:new int[]{1,2,8,64}){
				assertSame(expected,CladeBuilder.fromReads(reads,"fixture",threads));
				assertSame(expected,CladeBuilder.fromReads(new ListNum<Read>(reads,700),"fixture",threads));
			}
			final ArrayList<Read> empty=new ArrayList<Read>();assertSame(serial(empty,"empty"),CladeBuilder.fromReads(empty,"empty",4));
			final ArrayList<Read> zero=new ArrayList<Read>();zero.add(new Read(new byte[0],null,"zero",0));zero.add(reads.get(0));
			if(!entropy){assertSame(serial(zero,"zero"),CladeBuilder.fromReads(zero,"zero",2));}
			else{
				boolean serialRejected=false,parallelRejected=false;
				try{serial(zero,"zero");}catch(AssertionError e){serialRejected=true;}
				try{CladeBuilder.fromReads(zero,"zero",2);}catch(IllegalStateException e){parallelRejected=true;}
				check(serialRejected && parallelRejected,"EntropyTracker rejects zero-length sequences; factory must propagate the same failure rather than return partial counts");
			}
		}
		for(int i=0;i<reads.size();i++){check(Arrays.equals(before[i],reads.get(i).bases) && reads.get(i).id.equals(i<48?"contig"+i:i==48?"allN":"short"),"Construction must not mutate caller-owned Read objects");}
		CladeObject.calcCladeEntropy=false;
		final java.lang.management.ThreadMXBean threadBean=java.lang.management.ManagementFactory.getThreadMXBean();
		CladeBuilder.fromReads(reads,"warm-inline",1);
		final long started=threadBean.getTotalStartedThreadCount();
		CladeBuilder.fromReads(reads,"inline",1);
		check(threadBean.getTotalStartedThreadCount()==started,"t1 Read-list construction must spawn no thread");
		final Path input=root.resolve("fixture.fa");
		Files.write(input,">one\nACGTACGTACGTACGTACGTACGTACGT\n>two\nNNNNACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII));
		assertSame(CladeBuilder.fromFile(input.toString(),1),CladeBuilder.fromFile(input.toString(),4));
		boolean failed=false;try{CladeBuilder.fromReads(reads,"bad",0);}catch(IllegalArgumentException e){failed=true;}
		check(failed,"Invalid thread count must fail rather than silently skip work");
		System.out.println("PASS serial bytes, integer spectra, finished frequencies, sketch bytes, original-order entropy, Read/ListNum/file, t1/2/8/64, empty/N-rich/zero-length and immutable inputs");
	}
	static Clade serial(ArrayList<Read> reads,String name){
		final Clade c=new Clade(0,0,name);
		final EntropyTracker e=CladeObject.calcCladeEntropy?new EntropyTracker(CladeObject.entropyK,CladeObject.entropyWindow,false):null;
		for(Read r:reads){c.add(r.bases,e);}c.finish();return c;
	}
	static void assertSame(Clade a,Clade b){
		check(a.toBytes(new ByteBuilder()).toString().equals(b.toBytes(new ByteBuilder()).toString()),"Serialized Clade must equal unchanged serial construction");
		check(a.bases==b.bases && a.contigs==b.contigs && Float.floatToIntBits(a.entropy)==Float.floatToIntBits(b.entropy),"Base/contig totals and entropy bits must agree");
		for(int k=1;k<a.counts.length;k++){check(Arrays.equals(a.counts[k],b.counts[k]),"All integer spectra must agree at k="+k);}
		for(int k=0;k<a.frequencies.length;k++){check(Arrays.equals(a.frequencies[k],b.frequencies[k]),"All finished frequency arrays must agree at k="+k);}
		check((a.ddl==null)==(b.ddl==null),"DDL presence must be preserved");
		if(a.ddl!=null){check(a.ddl.toBytes(new ByteBuilder()).toString().equals(b.ddl.toBytes(new ByteBuilder()).toString()),"Serialized DDL buckets must agree");}
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
