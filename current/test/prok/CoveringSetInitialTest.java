package test.prok;

import prok.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import fileIO.ByteFile;
import map.LongHashSet;
import parse.LineParser1;
import shared.Shared;
import shared.Timer;

/** Fixed initial selection, caller-key projection and fixed-remainder stopping.
 * @author Raiden
 */
public final class CoveringSetInitialTest {
	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("cover-initial-");final Path pool=dir.resolve("pool.fa"),seed=dir.resolve("seed.fa");
		final StringBuilder p=new StringBuilder();for(int i=0;i<100;i++){p.append('>').append(i).append('\n').append(i<80?"AACC":i<90?"AACG":"TTTA").append('\n');}
		put(pool,p.toString());put(seed,">initial\nAACC\n>duplicate\nAACC\n");
		final int oldThreads=Shared.threads();
		try{
			final String a=run(dir,"one",pool,seed,true,1),b=run(dir,"four",pool,seed,true,4);
			check(Arrays.equals(Files.readAllBytes(dir.resolve("one.fa")),Files.readAllBytes(dir.resolve("four.fa"))),"Preloaded selection must be deterministic across parallel counting");
			check(a.equals(b),"Thread count must preserve all accounting");
			check(a.equals("100\t1\t2\t1\t3\t3\t2\t20\t10\ttrue\t1"),"One added design key must cover half of the fixed20-member remainder using caller3-mers: "+a);
			final String whole=run(dir,"whole",pool,seed,false,1);
			check(whole.equals("100\t1\t1\t0\t2\t3\t2\t20\t20\tfalse\t0"),"Whole-pool target and remainder target must remain distinct");
			final String design=new String(Files.readAllBytes(dir.resolve("one.design.fa")),StandardCharsets.US_ASCII);
			check(design.contains("\nAACC\n") && design.contains("\nAACG\n"),"Fixed original design key must survive augmentation unchanged");
			put(pool,">covered\nAACC\n");put(seed,">initial\nAACC\n>absent\nGGGG\n");
			final String complete=run(dir,"complete",pool,seed,true,1);
			check(complete.equals("1\t2\t2\t0\t3\t3\t2\t0\t0\ttrue\t0"),"Zero remainder must terminate immediately and retain keys absent from the pool");
			put(pool,">repeated\nAAAAAA\n");put(seed,">one\nAAAA\n");
			final String repeated=run(dir,"repeat",pool,seed,true,1);
			check(repeated.equals("1\t1\t1\t0\t1\t3\t2\t1\t1\ttrue\t0"),"Repeated occurrences and duplicate constituent keys must never become two distinct hits");
			check(CoveringSet.remainderTargetMet(1,100,.99) && !CoveringSet.remainderTargetMet(2,100,.99),"99% must apply to fixed initial remainder without floating-point boundary loss");
			check(!CoveringSet.remainderTargetMet(1,99,.99),"Small remainder must round required coverage upward");
			final LongHashSet forward=new LongHashSet(8),canonical=new LongHashSet(8);
			CoveringSet.addProjected(1L,4,3,false,forward);CoveringSet.addProjected(1L,4,3,true,canonical);
			check(forward.size()==2 && forward.contains(0) && forward.contains(1),"AAAC projects to AAA and AAC");
			check(canonical.size()==2 && canonical.contains(63) && canonical.contains(47),"Canonical projection must canonicalize each constituent independently");
		}finally{Shared.setThreads(oldThreads);}
		System.out.println("PASS CoveringSetInitialTest remainder projection distinct_hits fixed_keys threads empty_remainder");
	}
	private static String run(Path dir,String name,Path pool,Path seed,boolean remainder,int threads){
		new CoveringSet(new String[]{"in="+pool,"seedset="+seed,"out="+dir.resolve(name+".fa"),"outdesign="+dir.resolve(name+".design.fa"),"summary="+dir.resolve(name+".tsv"),"k=3","kdesign=4","coveragek=3","minhits=2","target=.5","targetremaining="+remainder,"steps=1","rcomp=f","t="+threads,"ow=f"}).process(new Timer());
		final ByteFile in=RrnaSourceTargetPrep.open(dir.resolve(name+".tsv").toString());
		try{check(in.nextLine()!=null,"Summary needs header");final byte[] row=in.nextLine();check(row!=null && in.nextLine()==null,"Exactly one summary row is required");return new String(row,StandardCharsets.US_ASCII);}finally{check(!in.close(),"Summary read failed");}
	}
	private static void put(Path path,String text) throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII));}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
