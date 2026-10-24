package test.prok;

import java.io.File;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import prok.NcrnaMitoRrnaExtract;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** End-to-end conservation and explicit expected sequences, including the circular join trap.
 * @author Ganyu */
public final class NcrnaMitoRrnaExtractTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Require fresh outdir=PATH");
		final String root=args[0].substring(7);require(new File(root).mkdir(),"Fixture directory already exists");
		final String roster="accession\ttaxid\torganelle_class\torganism\nAA_1.1\t11\tmitochondrion\tCafé α 🌿\nBB_2.1\t22\tkinetoplast\tfixture two\n";
		final ByteBuilder features=new ByteBuilder();
		feature(features,"2..5","12S α-rRNA café 🌿\t50%",false);
		feature(features,"complement(2..5)","16S ribosomal RNA",false);
		features.append("     rRNA            join(11..12,\n                     1..3)\n                     /product=\"café\n                     β 🌿\"\n".getBytes(java.nio.charset.StandardCharsets.UTF_8));
		feature(features,"complement(join(11..12,1..3))","18S ribosomal RNA",false);
		feature(features,"<2..>5","partial pseudo RNA",true);
		feature(features,"OTHER.1:1..4","remote",false);
		feature(features,"order(1..2,5..6)","ordered",false);
		feature(features,"4^5","between bases",false);
		feature(features,"6","single base",false);
		feature(features,"join(1..2,5..6)","gapped join",false);
		feature(features,"join(complement(11..12),complement(1..3))","individual complements",false);
		feature(features,"1..(3.4)","fuzzy right bound",false);
		final String first=record("AA_1.1",new String(features.array,0,features.length(),java.nio.charset.StandardCharsets.UTF_8)),second=record("BB_2.1","");
		write(root+"/roster.tsv",roster);write(root+"/in.gbff",first+second);
		run(root+"/in.gbff",root+"/roster.tsv",root+"/ok",2);
		final ObjectMap<String,String> candidates=fasta(root+"/ok.candidates.fa.gz"),genomes=fasta(root+"/ok.genomes.fa.gz");
		require(candidates.size()==8 && genomes.size()==2,"All supported rRNAs and every genome, including no-rRNA record, must survive");
		sequence(candidates,"AA_1.1_rrna_1","ACCG");sequence(candidates,"AA_1.1_rrna_2","CGGT");
		sequence(candidates,"AA_1.1_rrna_3","TCAAC");sequence(candidates,"AA_1.1_rrna_4","GTTGA");
		sequence(candidates,"AA_1.1_rrna_5","ACCG");sequence(candidates,"AA_1.1_rrna_9","G");
		sequence(candidates,"AA_1.1_rrna_10","AAGG");sequence(candidates,"AA_1.1_rrna_11","GAGTT");
		sequence(genomes,"BB_2.1","AACCGGTTAGTC");
		final String ledger=read(root+"/ok.loci.tsv.gz"),records=read(root+"/ok.records.tsv");
		require(lines(ledger)==13,"Every one of12 annotations, including unsupported ones, must have one ledger row");
		require(ledger.contains("\tUNSUPPORTED_REMOTE\t") && ledger.contains("\tUNSUPPORTED_ORDER\t") && ledger.contains("\tUNSUPPORTED_FUZZY\t"),"Unsupported locations cannot disappear or produce fabricated candidates");
		require(ledger.contains("\t<2..>5\t\"partial pseudo RNA\"\tEXTRACTED\ttrue\ttrue\t\"unitary\""),"Partial markers and pseudogene metadata must be retained");
		require(ledger.contains("\tjoin(11..12,1..3)\t\"café β 🌿\"\t") && ledger.contains("\t11..12:+,1..3:+\ttrue\t5\t"),"Wrapped UTF-8 product, raw location and circular segment order must survive");
		require(ledger.contains("\"12S α-rRNA café 🌿%0950%25\"") && records.contains("\tCafé α 🌿\t"),"UTF-8 metadata bytes must survive; TSV controls and percent must be escaped reversibly");
		require(records.contains("BB_2.1\t22\tkinetoplast\tfixture two\tNO_RRNA\t12\tcircular\t0\t0\t0\t"),"No-rRNA roster members remain available for CM fallback");
		require(ledger.contains("\t1..(3.4)\t\"fuzzy right bound\"\tUNSUPPORTED_FUZZY\t"),"Fuzzy range endpoints must be ledgered rather than silently lost or mistaken for malformed numeric coordinates");
		require(read(root+"/ok.summary.tsv").contains("\n2\t2\t2\t12\t8\t4\t1\t1\t1\t24\t0\n"),"Exact record/locus/status denominators must agree");
		write(root+"/duplicate.gbff",first+first);reject(root,"duplicate",root+"/duplicate.gbff",root+"/roster.tsv",2,"Repeated GBFF accession");
		write(root+"/missing.gbff",first);reject(root,"missing",root+"/missing.gbff",root+"/roster.tsv",2,"GBFF/roster conservation failed");
		require(read(root+"/missing.records.tsv").contains("\tMISSING\t"),"Missing records must be visible even on failure");
		write(root+"/duplicate_roster.tsv",roster+"AA_1.1\t11\tmitochondrion\tduplicate\n");
		reject(root,"duplicate_roster",root+"/in.gbff",root+"/duplicate_roster.tsv",3,"Duplicate roster accession");
		write(root+"/unknown.gbff",first+record("CC_3.1",""));reject(root,"unknown",root+"/unknown.gbff",root+"/roster.tsv",2,"GBFF accession absent from roster");
		write(root+"/length.gbff",(first+second).replace("12 bp","13 bp"));reject(root,"length",root+"/length.gbff",root+"/roster.tsv",2,"ORIGIN/LOCUS length mismatch");
		write(root+"/bounds.gbff",first.replace("2..5","2..99")+second);reject(root,"bounds",root+"/bounds.gbff",root+"/roster.tsv",2,"Location outside source sequence");
		write(root+"/truncated.gbff",first+second.substring(0,second.length()-3));reject(root,"truncated",root+"/truncated.gbff",root+"/roster.tsv",2,"Truncated GBFF record");
		final String fallbackFna=">BB_2.1 reordered full genome\nAACCGGTTAGTC\n>AA_1.1 reordered full genome\nAACCGGTTAGTC\n";
		final String contigFirst=contig(first),contigSecond=contig(second);
		write(root+"/fallback.fa",fallbackFna);write(root+"/contig.gbff",contigFirst+contigSecond);
		run(root+"/contig.gbff",root+"/roster.tsv",root+"/fallback",2,root+"/fallback.fa");
		final ObjectMap<String,String> fallback=fasta(root+"/fallback.candidates.fa.gz");
		require(fallback.size()==8,"CONTIG-only record with rRNAs must retain all supported candidates");
		for(String id:candidates.keys()){if(id!=null){sequence(fallback,id,candidates.get(id));}}
		require(fasta(root+"/fallback.genomes.fa.gz").size()==2,"CONTIG-only record without rRNA must still emit its full genome");
		require(read(root+"/fallback.summary.tsv").contains("\n2\t2\t2\t12\t8\t4\t1\t1\t1\t24\t2\n"),"Two accession-bound FASTA fallbacks must be counted exactly once");
		require(read(root+"/fallback.records.tsv").contains("\tFNA_CONTIG\t"+root+"/fallback.fa\tjoin(CHILD.1:1..12)\n"),"Fallback sequence source/path and raw CONTIG location must be retained");
		require(read(root+"/fallback.loci.tsv.gz").contains("\t1..(3.4)\t\"fuzzy right bound\"\tUNSUPPORTED_FUZZY\t"),"CONTIG must finish the last feature; its assembly cannot become a product continuation");
		write(root+"/mixed.gbff",first+contigSecond);
		run(root+"/mixed.gbff",root+"/roster.tsv",root+"/mixed",2,root+"/fallback.fa");
		require(read(root+"/mixed.summary.tsv").contains("\n2\t2\t2\t12\t8\t4\t1\t1\t1\t24\t1\n") && read(root+"/mixed.records.tsv").contains("\tGBFF_ORIGIN\t"),"Mixed ORIGIN/CONTIG input cannot duplicate or skip either source");
		reject(root,"no_fna",root+"/contig.gbff",root+"/roster.tsv",2,"CONTIG record requires fna=");
		write(root+"/missing_fna.fa",">BB_2.1\nAACCGGTTAGTC\n");
		reject(root,"missing_fna",root+"/contig.gbff",root+"/roster.tsv",2,"FASTA/roster conservation failed",root+"/missing_fna.fa");
		write(root+"/duplicate_fna.fa",fallbackFna+">AA_1.1\nAACCGGTTAGTC\n");
		reject(root,"duplicate_fna",root+"/contig.gbff",root+"/roster.tsv",2,"Repeated FASTA accession",root+"/duplicate_fna.fa");
		write(root+"/wrong_length_fna.fa",fallbackFna.replace("AACCGGTTAGTC","AACCGGTTAGT"));
		reject(root,"wrong_length_fna",root+"/contig.gbff",root+"/roster.tsv",2,"FASTA/LOCUS length mismatch",root+"/wrong_length_fna.fa");
		write(root+"/unknown_fna.fa",fallbackFna.replace("BB_2.1","CC_3.1"));
		reject(root,"unknown_fna",root+"/contig.gbff",root+"/roster.tsv",2,"FASTA accession absent from roster",root+"/unknown_fna.fa");
		System.out.println("MITO_RRNA_EXTRACT_TEST_PASS records=2 features=12 extracted=8 unsupported=4 negative=12 circular_join=true partial_pseudo=true no_rRNA_retained=true contig_fallback=true");
	}
	private static void reject(String root,String name,String gbff,String roster,int expected,String message){
		reject(root,name,gbff,roster,expected,message,null);
	}
	private static void reject(String root,String name,String gbff,String roster,int expected,String message,String fna){
		boolean rejected=false;try{run(gbff,roster,root+"/"+name,expected,fna);}catch(IllegalArgumentException e){require(e.getMessage().contains(message),"Unexpected failure for "+name+": "+e);rejected=true;}
		require(rejected,"Expected fail-loud fixture: "+name);require(!new File(root+"/"+name+".summary.tsv").exists(),"A failed census cannot leave a success summary");
	}
	private static void run(String gbff,String roster,String prefix,int expected){run(gbff,roster,prefix,expected,null);}
	private static void run(String gbff,String roster,String prefix,int expected,String fna){
		String[] args={"in="+gbff,"roster="+roster,"outprefix="+prefix,"expectedrecords="+expected,"t=1","zl=9"};
		if(fna!=null){args=java.util.Arrays.copyOf(args,7);args[6]="fna="+fna;}NcrnaMitoRrnaExtract.main(args);
	}
	private static String contig(String record){return record.replace("ORIGIN\n        1 aaccggttag tc\n","CONTIG      join(CHILD.1:1..12)\n");}
	private static String record(String accession,String features){return "LOCUS       "+accession+" 12 bp DNA circular\nACCESSION   "+accession.substring(0,accession.indexOf('.'))+"\nVERSION     "+accession+"\nFEATURES             Location/Qualifiers\n     source          1..12\n"+features+"ORIGIN\n        1 aaccggttag tc\n//\n";}
	private static void feature(ByteBuilder b,String location,String product,boolean pseudo){b.append("     rRNA            ").append(location).append("\n                     /product=\"").append(product.getBytes(java.nio.charset.StandardCharsets.UTF_8)).append("\"\n");if(pseudo){b.append("                     /pseudo\n                     /pseudogene=\"unitary\"\n");}}
	private static void sequence(ObjectMap<String,String> map,String id,String expected){require(expected.equals(map.get(id)),"Wrong extracted bases for "+id+": "+map.get(id)+" expected="+expected);}
	private static int lines(String s){int n=0;for(int i=0;i<s.length();i++){if(s.charAt(i)=='\n'){n++;}}return n;}
	private static ObjectMap<String,String> fasta(String path){
		final ObjectMap<String,String> result=new ObjectMap<String,String>(String.class,String.class);
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		st.start();try{for(ListNum<Read> batch;(batch=st.nextList())!=null;){for(Read r:batch.list){
			final String id=r.id.split(" ",2)[0];require(!result.contains(id),"Duplicate FASTA identifier");
			result.put(id,new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII));
		}st.returnList(batch);}require(!st.errorState(),"Fixture FASTA read error");}finally{require(!ReadWrite.closeStream(st),"Fixture FASTA close error");}return result;
	}
	private static String read(String path){final ByteFile f=ByteFile.makeByteFile(path,true);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=f.nextLine())!=null;){b.append(line).nl();}}finally{require(!f.close(),"Fixture read failure");}return new String(b.array,0,b.length(),java.nio.charset.StandardCharsets.UTF_8);}
	private static void write(String path,String text){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();w.print(text.getBytes(java.nio.charset.StandardCharsets.UTF_8));require(!w.poisonAndWait(),"Fixture write failure");}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
