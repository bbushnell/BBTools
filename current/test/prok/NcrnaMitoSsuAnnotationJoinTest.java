package test.prok;

import java.io.File;
import java.util.Arrays;
import prok.NcrnaMitoSsuAnnotationJoin;

/** Position-preserving annotation projection and independent name/CM evidence.
 * @author Brian, Ganyu */
public final class NcrnaMitoSsuAnnotationJoinTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be fresh");
		final String[] ch=new String[32];Arrays.fill(ch,"unused");ch[0]="case_id";ch[10]="any_supported_bp";ch[31]="alignment_identity_min";
		write(root,"calls.tsv",String.join("\t",ch)+"\n"+call("F",101,"+",10)+call("R",201,"+",10)+call("J",301,"+",0)
			+call("Z",401,"+",0)+call("M",501,"+",0)+call("T",601,"-",10)+call("N",701,"+",10));
		final String[] ah=new String[28];Arrays.fill(ah,"unused");ah[0]="locus_id";ah[8]="status";ah[23]="source_segments_in_transcript_order";ah[26]="source_gbff";
		write(root,"loci.tsv",String.join("\t",ah)+"\n"
			+annotation("SSU","RF","SSU_SUPPORT_REVIEW","\"12S ribosomal RNA\"","+",11,23,"11..12:+,22..23:+",4)
			+annotation("GAP","RF","SSU_SUPPORT_REVIEW","small subunit ribosomal RNA","+",15,19,"15..19:+",5)
			+annotation("OPPOSITE","RF","SSU_SUPPORT_REVIEW","12S ribosomal RNA","-",10,25,"10..25:-",16)
			+annotation("CONFLICT","RF","SSU_LSU_CONFLICT_REVIEW","small-subunit ribosomal RNA","+",20,21,"20..21:+",2)
			+annotation("LSU","RF","LSU_SUPPORT_REVIEW","16S ribosomal RNA","+",24,25,"24..25:+",2)
			+annotation("UNRESOLVED","RF","NO_SAME_STRAND_GA_HIT_UNRESOLVED","16S ribosomal RNA","+",10,10,"10..10:+",1)
			+annotation("REVERSE","RR","SSU_SUPPORT_REVIEW","12S ribosomal RNA","-",41,42,"41..42:-",2)
			+annotation("WRONG_STRAND","RR","SSU_SUPPORT_REVIEW","12S ribosomal RNA","+",40,44,"40..44:+",5));
		final String[] hh=new String[19];Arrays.fill(hh,"unused");hh[0]="case_id";hh[14]="original_supported_intervals0";hh[16]="reference_supported_intervals1";hh[18]="passes_alignment_identity_076";final String header=String.join("\t",hh)+"\n";
		final String forward=alignment("F","RF",0,"100-103,104-109","10-13,20-25",10,0,true);
		final String rest=alignment("F","RF",256,"100-103,104-109","10-13,20-25",10,0,true)
			+alignment("R","RR",16,"205-209,200-204","40-44,50-54",10,0,true)
			+alignment("J","RF",0,"300-309","10-19",10,0,false)
			+alignment("M","RF",0,"","",0,10,true)
			+alignment("T","RF",16,"605-609,600-604","10-14,20-24",10,0,true)
			+alignment("N","NO_ANNOTATION",0,"700-709","10-19",10,0,true);
		write(root,"alignments.tsv",header+forward+rest);final String pins=NcrnaMitoCrossLibraryTest.bind(root,sha,"calls.tsv","loci.tsv","alignments.tsv");run(root,pins,"good","alignments.tsv");
		final String calls=read(root+"/good.calls.tsv");check(calls.split("\n").length==8,"All seven original calls must remain, including zero-hit and unnamed-reference calls");
		final String[] f=NcrnaMitoCrossLibraryTest.find(calls,"F\t"),r=NcrnaMitoCrossLibraryTest.find(calls,"R\t"),t=NcrnaMitoCrossLibraryTest.find(calls,"T\t");
		checkSupport(f,0,4);checkSupport(f,1,2);checkSupport(f,2,2);checkSupport(f,3,1);checkSupport(f,4,6);checkSupport(f,5,2);
		checkSupport(r,0,2);checkSupport(r,4,2);checkSupport(t,0,4);checkSupport(t,4,6);
		for(String id:new String[]{"J","Z","M","N"}){final String[] row=NcrnaMitoCrossLibraryTest.find(calls,id+'\t');for(int i=0;i<6;i++){checkSupport(row,i,0);}}
		final String annotations=read(root+"/good.annotations.tsv.gz");
		check(!annotations.contains("\tGAP\t"),"Neither annotation holes nor reference deletion gaps supply paired query support");
		final String[] fs=NcrnaMitoCrossLibraryTest.find(annotations,"F\tCF\tRF\tSSU\t"),rs=NcrnaMitoCrossLibraryTest.find(annotations,"R\tCR\tRR\tREVERSE\t");
		check(fs[11].equals("4") && fs[14].equals("2") && fs[18].equals("101-102,101-102,106-107,106-107"),"Secondary duplicate overlaps remain auditable but unique query coverage stays four bases");
		check(rs[11].equals("2") && rs[18].equals("207-208"),"Reverse annotation intersections must project back from the high query endpoint");
		final String[] opposite=NcrnaMitoCrossLibraryTest.find(annotations,"F\tCF\tRF\tOPPOSITE\t");check(opposite[11].equals("0") && opposite[13].equals("10"),"Wrong-strand geometric overlap must remain separate from support");
		final String[] lsu=NcrnaMitoCrossLibraryTest.find(annotations,"F\tCF\tRF\tLSU\t");check(lsu[9].equals("CM_LSU_ONLY") && lsu[10].equals("no"),"Mitochondrial 16S naming cannot silently become SSU support");
		check(NcrnaMitoSsuAnnotationJoin.namedSsu("\"12S ribosomal RNA\"") && NcrnaMitoSsuAnnotationJoin.namedSsu("small subunit ribosomal RNA")
			&& !NcrnaMitoSsuAnnotationJoin.namedSsu("16S ribosomal RNA") && !NcrnaMitoSsuAnnotationJoin.namedSsu("predicted 12S ribosomal RNA"),"Name-only evidence follows exact12S or explicit small-subunit wording");
		write(root,"bad.tsv",header+alignment("F","RF",0,"100-103,104-109","10-12,20-25",10,0,true)+rest);
		final String badPins=NcrnaMitoCrossLibraryTest.bind(root,sha,"calls.tsv","loci.tsv","bad.tsv");boolean failed=false;
		try{run(root,badPins,"bad","bad.tsv");}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+"/bad.calls.tsv").exists(),"Unequal query/reference pair lengths must fail before result emission");
		System.out.println("MITO_SSU_ANNOTATION_JOIN_TEST_PASS forward reverse call_orientation source_segments deletion_gaps duplicate_union rejected zero_hit ambiguous16S independent_name_CM invalid_pair");
	}
	private static void run(String root,String pins,String out,String alignments){NcrnaMitoSsuAnnotationJoin.main(new String[]{"calls="+root+"/calls.tsv","loci="+root+"/loci.tsv","alignments="+root+"/"+alignments,"verified="+pins,"outprefix="+root+"/"+out,"expected=7"});}
	private static void checkSupport(String[] row,int category,int count){final int index=32+3*category;check(row[index].equals(count>0?"yes":"no") && Integer.parseInt(row[index+1])==count && Math.abs(Double.parseDouble(row[index+2])-count/10.0)<1e-9,"Source-annotation support must count unique known query bases: "+row[0]+" category="+category);}
	private static String call(String id,int start,String strand,int known){final String[] p=new String[32];Arrays.fill(p,"NA");p[0]=id;p[1]="C"+id;p[5]=""+start;p[6]=""+(start+9);p[7]=strand;p[8]="10";p[10]=""+known;p[31]="0.76";return String.join("\t",p)+"\n";}
	private static String alignment(String id,String ref,int flag,String query,String reference,int matches,int unknown,boolean passes){final String[] p=new String[19];Arrays.fill(p,"NA");p[0]=id;p[1]="C"+id;p[2]="tile";p[3]=ref;p[4]=""+flag;p[5]=""+((flag&256)==0);p[6]=""+matches;p[7]="0";p[11]=""+unknown;p[14]=query;p[16]=reference;p[18]=""+passes;return String.join("\t",p)+"\n";}
	private static String annotation(String id,String ref,String status,String product,String strand,int min,int max,String segments,int length){final String[] p=new String[28];Arrays.fill(p,"NA");p[0]=id;p[1]=ref;p[8]=status;p[9]=id;p[10]=ref;p[14]=segments;p[15]=product;p[16]="EXTRACTED";p[17]="false";p[18]="false";p[20]=strand;p[21]=""+min;p[22]=""+max;p[23]=segments;p[25]=""+length;p[26]="fixture.gbff";return String.join("\t",p)+"\n";}
	private static void write(String root,String name,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+name,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
