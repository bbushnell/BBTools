package test.prok;

import prok.*;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import parse.LineParser1;
import shared.Shared;
import stream.Read;

/** Three/four distinct keys, ambiguity and repeated copies must remain distinct concepts.
 * @author Raiden
 */
public final class Euk18sSeedMinimumCoverageTest {
	public static void main(String[] args){
		check(args.length==1 && new File(args[0]).mkdirs(),"Fresh fixture directory required");
		final String dir=args[0];Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		final LongHashSet seeds=new LongHashSet();seeds.add(0);seeds.add(1);seeds.add(2);seeds.add(3);
		final Euk18sSeedMinimumCoverage tool=new Euk18sSeedMinimumCoverage(seeds);
		final String a="AAAAAAAAAAAAAAAAA",c="AAAAAAAAAAAAAAAAC",g="AAAAAAAAAAAAAAAAG",t="AAAAAAAAAAAAAAAAT";
		final String three=a+"N"+c+"N"+g,four=three+"N"+t;
		check(tool.count(bytes(a+a+a))==1,"Repeated occurrences cannot satisfy a distinct-key minimum");
		check(tool.count(bytes(three))==3 && tool.count(bytes(four))==4,"Independent forward keys must reach their exact minima");
		check(tool.count(bytes("AAAAAAAANAAAAAAAAA"))==0,"Ambiguity breaks keys; flanks cannot concatenate");
		check(tool.count(bytes("TTTTTTTTTTTTTTTTT"))==0,"Forward-only coverage must not introduce reverse-complement matches");
		check(tool.count(bytes(a+"CG"))==2,"Overlapping keys remain distinct while absent overlaps do not count");
		final ByteStreamWriter out=RrnaSourceTargetPrep.writer(dir+"/members.fa");
		try{out.println(">same assembly=A\n"+three+"\n>same assembly=B\n"+three+"\n>four assembly=C\n"+four+"\n>repeat assembly=A\n"+a+a);}
		finally{RrnaCmWindowParity.finish(out);}
		tool.run(dir+"/members.fa","fixture",dir+"/min3",3,4);
		tool.run(dir+"/members.fa","fixture",dir+"/min4",4,4);
		checkSummary(dir+"/min3.summary.tsv",3,3,1);checkSummary(dir+"/min4.summary.tsv",4,1,3);
		final ByteFile below=RrnaSourceTargetPrep.open(dir+"/min4.below.tsv");final LineParser1 p=new LineParser1('\t');
		try{check(below.nextLine()!=null,"Below-min header required");p.set(below.nextLine());check(p.termEquals("A",2)&&p.parseInt(6)==3,"First three-key member must stay below four");p.set(below.nextLine());check(p.termEquals("B",2)&&p.parseInt(6)==3,"Identical sequence and name from another assembly must remain another copy");p.set(below.nextLine());check(p.parseInt(6)==1,"Repeated single key must remain one");check(below.nextLine()==null,"Exactly three below-min members expected");}finally{check(!below.close(),"Below-min fixture read failed");}
		final int[] lengths={915,916,1599,1600,2200,2201,2500,2501},classes={1,2,2,3,3,4,4,5};
		for(int i=0;i<lengths.length;i++){check(Euk18sSeedMinimumCoverage.lengthClass(lengths[i])==classes[i],"Length-bin boundary changed");}
		boolean failed=false;try{tool.run(dir+"/members.fa","fixture",dir+"/wrong",3,5);}catch(IllegalArgumentException e){failed=true;}
		check(failed,"Wrong corpus denominator must fail rather than publish coverage");
		System.out.println("EUK18S_SEED_MINIMUM_TEST_PASS distinct3_4 repeats ambiguity orientation raw_copy_conservation length_bins denominator");
	}
	static void checkSummary(String path,int minimum,int covered,int missing){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{check(in.nextLine()!=null,"Summary header required");p.set(in.nextLine());check(p.termEquals("ALL",1)&&p.parseInt(2)==minimum&&p.parseInt(4)==4&&p.parseInt(5)==covered&&p.parseInt(6)==missing,"Coverage must partition all four raw members at the declared minimum");}finally{check(!in.close(),"Summary read failed");}
	}
	static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
