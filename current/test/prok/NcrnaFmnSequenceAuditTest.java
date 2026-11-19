package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaConstructionMembers;
import prok.NcrnaFmnSequenceAudit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal extraction, input binding and independent exhaustive substring oracle.
 * @author Ganyu
 */
public final class NcrnaFmnSequenceAuditTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFmnSequenceAuditTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown sequence-audit test option");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Use out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture output must be fresh");
			matcher();
			fixture(root.getPath());
			pp.outstream.println("NCRNA_FMN_SEQUENCE_AUDIT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private static void matcher(){
		final byte[][] queries={bytes(Q), bytes(P)};
		final NcrnaFmnSequenceAudit.Matcher matcher=new NcrnaFmnSequenceAudit.Matcher(queries);
		final byte[] reverse=bytes(P);
		AminoAcid.reverseComplementBasesInPlace(reverse);
		final byte[][] targets={bytes("NN"+Q+"NN"), reverse, bytes(Q.toLowerCase(Locale.ROOT)), bytes(Q.substring(0, 30)),
			bytes(Q.substring(0, 20)+"NN"+Q.substring(20)), bytes("NNNNNN"), bytes(P+Q), bytes("")};
		for(final byte[] target:targets){
			final NcrnaFmnSequenceAudit.Match observed=new NcrnaFmnSequenceAudit.Match();
			matcher.scan(target, "fixture", observed);
			int expected=0;
			for(final byte[] q:queries){
				final byte[] rc=q.clone();
				AminoAcid.reverseComplementBasesInPlace(rc);
				expected=Math.max(expected, Math.max(longest(target, q), longest(target, rc)));
			}
			if(expected<31){expected=0;}
			check(observed.length==expected, "Exact seed extension must equal an independent all-position substring oracle");
		}
		final NcrnaFmnSequenceAudit.Match split=new NcrnaFmnSequenceAudit.Match();
		matcher.scan(bytes(Q.substring(0, 25)), "left", split);
		matcher.scan(bytes(Q.substring(25)), "right", split);
		check(split.length==0, "Exact seeds may not bridge separate contigs");
		final NcrnaFmnSequenceAudit.Match minus=new NcrnaFmnSequenceAudit.Match();
		matcher.scan(reverse, "minus", minus);
		check(minus.length==P.length() && minus.query==1 && minus.minus, "A reverse-strand full match must retain original query identity");
	}
	/** Independent exhaustive substring comparisons; no kmer packing or index. */
	private static int longest(final byte[] target, final byte[] query){
		int best=0;
		for(int a=0; a<target.length; a++){
			for(int b=0; b<query.length; b++){
				int n=0;
				while(a+n<target.length && b+n<query.length && upper(target[a+n])==upper(query[b+n])){n++;}
				best=Math.max(best, n);
			}
		}
		return best;
	}
	private static void fixture(final String root){
		final String[] paths={root+"/eval0.fa", root+"/eval1.fa", root+"/discovery.fa", root+"/elizabeth.fa", root+"/empedo.fa"};
		final byte[] reverse=bytes(P);
		AminoAcid.reverseComplementBasesInPlace(reverse);
		final String[] bases={"NN"+Q+"NN", "NN"+new String(reverse, StandardCharsets.US_ASCII)+"NN", Q.toLowerCase(Locale.ROOT), "NNNNNNNN", Q.substring(0, 20)};
		for(int i=0; i<paths.length; i++){write(paths[i], ">contig"+i+" description\n"+bases[i]+"\n");}
		String lineages="cohort\tinput_id\tfasta\tassembly_group\tassembly\tclass\torder\tfamily\tgenus\ttaxonomy_status\tCM_admitted_cores\n";
		String metadata="accession\tncbi_genome_category\tncbi_submitter\tncbi_assembly_level\tcheckm2_completeness\tgenome_size\tcontig_count\tncbi_total_length\tncbi_contig_count\tncbi_genome_representation\tgc_percentage\n";
		String census="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n";
		for(int i=0; i<5; i++){
			final String cohort=i<2 ? "evaluation" : "discovery";
			final String genus=i<3 ? "Chryseobacterium" : i==3 ? "Elizabethkingia" : "Empedobacter";
			lineages+=cohort+'\t'+(i+1)+'\t'+paths[i]+'\t'+i+"\tGCF_"+i+".1\tClass\tOrder\tFamily\t"+genus+"\tGTDB_ACCESSION_MATCH\t"+(i<2 ? 1 : 0)+'\n';
			metadata+="RS_GCF_"+i+".1\tnone\tLab α\tContig\t99\t"+bases[i].length()+"\t1\t"+bases[i].length()+"\t1\tfull\t40\n";
			if(i>=2){census+=(i+1)+"\t"+paths[i]+'\t'+new File(paths[i]).length()+'\t'+bases[i].length()+"\t1\n";}
		}
		String cores=NcrnaConstructionMembers.EXTRACT_HEADER+'\n';
		for(int i=0; i<2; i++){
			final int length=i==0 ? Q.length() : P.length();
			cores+="q"+i+"\tevaluation\t"+(i+1)+'\t'+paths[i]+'\t'+i+"\t1\tbacteria\tBacteroidota\tcontig"+i+'\t'+(i==0 ? "+" : "-")+"\t3-"+(length+2)+'\t'+length+"\tCM_FULL_CLEAN\tEVALUATION_ASSEMBLY_HELD_OUT\ttrue\ttrue\tfalse\tfalse\tNA\tevaluation:1\tPENDING_SEQUENCE_AND_FLANK_VALIDATION\n";
		}
		write(root+"/lineages.tsv", lineages);
		write(root+"/metadata.tsv", metadata);
		write(root+"/census.tsv", census);
		write(root+"/cores.tsv", cores);
		write(root+"/previous.tsv", NcrnaConstructionMembers.CROSS_HEADER+'\n'+"discovery\t3\t"+paths[2]+"\t2\tGCF_2.1\t1\tbacteria\tBacteroidota\tClass\tOrder\tGTDB_ACCESSION_MATCH\t0\n");
		final NcrnaFmnSequenceAudit.Config c=new NcrnaFmnSequenceAudit.Config();
		c.lineages=root+"/lineages.tsv";
		c.metadata=root+"/metadata.tsv";
		c.census=root+"/census.tsv";
		c.cores=root+"/cores.tsv";
		c.previous=root+"/previous.tsv";
		c.output=root+"/audit";
		c.queries=2;
		c.discovery=1;
		NcrnaFmnSequenceAudit.run(c);
		check(rows(c.output+".queries.fa").equals(java.util.Arrays.asList(">q0", Q, ">q1", P)), "Exact extraction must respect inclusive source coordinates and minus-strand orientation");
		check(rows(c.output+".summary.tsv").get(1).equals("2\t3\t1\t1\tEXACT_SEQUENCE_AUDIT_NOT_CM_OR_SPLIT"), "Summary must include zero-hit controls and actual previous-selection overlap");
		final String[] matched=rows(c.output+".presence.tsv").get(1).split("\t", -1);
		check(matched[7].equals(Integer.toString(Q.length())) && matched[13].equals(Integer.toString(Q.length())), "Actual lowercase count and best contiguous length must match the original input bytes");
		c.census=root+"/bad_census.tsv";
		write(c.census, census.replace("\t"+Q.length()+"\t1\n", "\t999\t1\n"));
		c.output=root+"/bad";
		boolean rejected=false;
		try{NcrnaFmnSequenceAudit.run(c);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("dimensions changed");}
		check(rejected, "Changed original scan dimensions must stop the diagnostic");
	}
	private static void write(final String path, final String value){try(Output out=new Output(path)){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}}
	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}}
		return result;
	}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static byte upper(final byte b){return b>='a' && b<='z' ? (byte)(b-32) : b;}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private static final String Q="ACGTTGCAACCTGATCGATGCTAGCTTACCGGATCGTACGA";
	private static final String P="GGCATACCTGATAGTCGATCGTTAACGGCATTCGATACCGTAGCTA";
}
