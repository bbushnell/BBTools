package test.prok;

import prok.*;

import java.io.File;
import java.util.ArrayList;
import fileIO.ByteStreamWriter;
import map.ObjectSet;
import shared.Timer;

/** Paper-checkable avoidance, coverage and projection fixtures using the actual selector.
 * @author Raiden */
public final class CoveringSetAvoidTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Fresh fixture outdir required");
		final String dir=args[0].substring(7);check(new File(dir).mkdir(),"Do not overwrite earlier test evidence");
		final String main=dir+"/pool.fa",aa=dir+"/aa.fa.gz",cc=dir+"/cc.fa",tie=dir+"/tie.fa",empty=dir+"/empty.fa";
		write(main,">a\nAA\n>b\nCC\n>c\nGG\n>d\nTT\n");write(aa,">a\nAAAA\n");write(cc,">c\nCCCC\n");
		write(tie,">a\nAA\n>c\nCC\n");write(empty,">n\nNNNNN\n");
		final String baseline=run(dir,"baseline",main),disabled=run(dir,"empty",main,"avoid="+empty,"avoidcount=1");
		check(read(baseline).equals(read(disabled)),"Zero valid avoidance keys must preserve exact baseline FASTA bytes");
		check(matches(run(dir,"union",main,"avoid="+aa+","+cc,"avoidcount=1"),"GG","TT"),"Each file contributes its own top1; union bans AA and CC, even when target becomes unreachable");
		check(matches(run(dir,"tie",main,"avoid="+tie,"avoidcount=1"),"CC","GG","TT"),"Equal occurrence counts choose the smaller encoded key AA deterministically");
		check(matches(run(dir,"canonical",main,"avoid="+aa,"avoidcount=1","rcomp=t"),"GG"),"Canonical avoidance must ban both AA and TT while retaining the GG/CC equivalence class");
		final String multi=dir+"/multihit.fa";write(multi,">a\nAACG\n>b\nAACT\n");
		check(matches(run(dir,"multihit",multi,"avoid="+aa,"avoidcount=1","minhits=2"),"AC","CG","CT"),"Banned AA must not enter selectedSet and falsely satisfy the two-hit eviction threshold");
		final String once=dir+"/once.fa";write(once,">a\nAA\n");
		check(matches(run(dir,"floor",multi,"avoid="+once,"avoidcount=1","mincount=2"),"AC"),"Main-pool mincount cannot suppress a low-count avoidance key before its independent ranking");
		final String design=dir+"/design.fa";write(design,">a\nAAACG\n");
		check(matches(run(dir,"design",design,"avoid="+aa,"avoidcount=1","kdesign=3"),"AC","CG"),"A design word containing banned output AA must be excluded before projection, not filtered only afterward");
		final String large=dir+"/parallel.fa";final StringBuilder seq=new StringBuilder();
		for(int i=0;i<80;i++){seq.append('>').append(i).append('\n').append((i&1)==0?"AACG":"AACT").append('\n');}write(large,seq.toString());
		final String one=run(dir,"one",large,"avoid="+aa,"avoidcount=1","minhits=2");
		final String many=run(dir,"many",large,"avoid="+aa,"avoidcount=1","minhits=2","t=8");
		check(read(one).equals(read(many)),"Partition count and worker count must not change selected key order");
		check(matches(many,"AC","CG","CT"),"Bounded top2 candidate heap must keep progressing after the high-count banned and already-selected keys");
		final String initial=dir+"/initial.fa";write(initial,">fixed\nAC\n");
		check(matches(run(dir,"augment",multi,"avoid="+aa,"avoidcount=1","minhits=2","seedset="+initial,"targetremaining=t"),"AC","CG","CT"),"Permitted fixed seeds survive augmentation with avoidance and two-hit coverage");
		boolean conflict=false;try{run(dir,"conflict",multi,"avoid="+aa,"avoidcount=1","seedset="+once);}catch(IllegalArgumentException e){conflict=e.getMessage().contains("conflicts with avoid");}
		check(conflict && !new File(dir+"/conflict.out.fa").exists(),"Conflicting fixed seeds and bans must fail before output, never silently drop a fixed seed");
		expectInvalid(main,dir+"/invalid.fa","avoidcount=1");
		expectInvalid(main,dir+"/invalid.fa","avoid="+aa,"avoidcount=0");
		expectInvalid(main,dir+"/invalid.fa","avoid="+aa,"alphabet=amino");
		System.out.println("PASS per-file top-N union, deterministic ties, compressed input, rcomp, multi-hit coverage, mincount independence, design projection, threaded parity, zero-key control and invalid options");
	}
	private static String run(String dir,String tag,String in,String... extra){
		final String out=dir+"/"+tag+".out.fa";final ArrayList<String> args=new ArrayList<String>();
		for(String a:new String[]{"in="+in,"out="+out,"k=2","step=1","target=1","rcomp=f","t=1","ow=f"}){args.add(a);}
		for(String a:extra){args.add(a);}new CoveringSet(args.toArray(new String[args.size()])).process(new Timer());return out;
	}
	private static ObjectSet<String> words(String file){
		final ObjectSet<String> out=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(file,r->{check(out.add(new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII)),"Selected output words must be unique");});return out;
	}
	private static boolean matches(String file,String... expected){final ObjectSet<String> actual=words(file);if(actual.size()!=expected.length){return false;}for(String word:expected){if(!actual.contains(word)){return false;}}return true;}
	private static String read(String path){
		final fileIO.ByteFile in=RrnaResourceIO.open(path);final StringBuilder out=new StringBuilder();
		try{for(byte[] line;(line=in.nextLine())!=null;){out.append(new String(line,java.nio.charset.StandardCharsets.US_ASCII)).append('\n');}}
		finally{check(!in.close(),"Fixture output read failed");}return out.toString();
	}
	private static void expectInvalid(String in,String out,String... options){
		final ArrayList<String> args=new ArrayList<String>();args.add("in="+in);args.add("out="+out);for(String x:options){args.add(x);}
		boolean failed=false;try{new CoveringSet(args.toArray(new String[args.size()]));}catch(IllegalArgumentException e){failed=true;}
		check(failed,"Invalid avoidance configuration must fail before any selection");
	}
	private static void write(String path,String text){final ByteStreamWriter out=RrnaResourceIO.writer(path);out.print(text);RrnaResourceIO.finish(out);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
