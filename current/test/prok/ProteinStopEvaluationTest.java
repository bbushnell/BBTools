package test.prok;

import prok.ProteinStopEvaluation;
import prok.ProteinSupportCompare;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;

import gff.GffLine;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;

/** Literal independent coordinate/strand/exception oracles; no model dependency.
 * @author Sayu
 */
public final class ProteinStopEvaluationTest {

	public static void main(String[] args){
		one("ATGAAATAA", "+", 0, ".", "primary");
		one("TTATTTCAT", "-", 0, ".", "primary");
		one("AAAAAATAA", "+", 0, ".", "primary"); // start policy is not stop eligibility
		one("GTGAAATAA", "+", 0, "partial=10", "primary");
		one("TTATTTCAC", "-", 0, "partial=01", "primary");
		one("AGTGAAATAA", "+", 1, "partial=10", "primary");
		one("TTATTTCACT", "-", 1, "partial=01", "primary");
		one("CCGTGAAATAA", "+", 2, "partial=10", "primary");
		one("TTATTTCACGG", "-", 2, "partial=01", "primary");
		one("ATGAAAA", "+", 0, "partial=01", "partial_3prime");
		one("TTTTCAT", "-", 0, "partial=10", "partial_3prime");
		one("ATGAAATAA", "+", 0, "partial=true", "partial_side_unknown");
		one("ATGAAATAA", "+", 0, "partial=true;start_range=.,1", "uncertain_range");
		one("ATGAAATAA", "+", -1, ".", "unknown_phase");
		one("ATGAAATAA", ".", 0, ".", "unknown_strand");
		one("ATGAAATAA", "+", 0, "transl_table=2", "unsupported_table");
		one("ATGTGATAA", "+", 0, "transl_table=4", "primary");
		one("ATGTGATAA", "+", 0, "transl_table=25", "primary");
		one("ATGTGATAA", "+", 0, ".", "internal_stop");
		one("ATGAAATGA", "+", 0, "transl_table=4", "missing_terminal_stop");
		one("ATGNNNTAA", "+", 0, ".", "ambiguous_codon");
		one("ATGAAATAA", "+", 0, "transl_except=(pos:4..6%2Caa:Sec)", "translation_exception");
		one("ATGAAATAA", "+", 0, "exception=ribosomal slippage", "translation_exception");
		one("ATGAAATAA", "+", 0, "pseudo=true", "pseudo");
		one("ATGAAATAA", "+", 0, "pseudogene=unitary", "pseudo");
		one("ATGAAATAA", "+", 0, "Parent=absent", "missing_parent");
		one("ATGAAATAA", "+", 0, "Parent=a,b", "multiple_parents");
		one("ATGAAATAA", "+", 0, "Note=pseudo=true;partial=00", "primary");
		duplicatesAndMasks(); parentAndMultipart(); circularOrigin(); nativeSynthetic(); failures();
		if(args.length>0){fileRoundTrip(args[0]);}
		System.err.println("STOP_EVALUATION_TEST_PASS literal_rows=28 duplicates_masks=true parents_multipart=true invalid_inputs=true");
	}

	private static void fileRoundTrip(String dir){
		final String fasta=dir+"/genome.fa", ref=dir+"/reference.gff", before=dir+"/before.gff", after=dir+"/after.gff";
		write(fasta, ">ctg\nATGATGAAATAAATGAAATAA");
		write(ref, "ctg\tfixture\tCDS\t1\t12\t.\t+\t0\tID=x\nctg\tfixture\tCDS\t1\t12\t.\t+\t0\tID=x");
		write(before, "ctg\tfixture\tCDS\t1\t12\t.\t+\t0\t.");
		write(after, "ctg\tfixture\tCDS\t1\t12\t.\t+\t0\t.\nctg\tfixture\tCDS\t4\t12\t.\t+\t0\t.\nctg\tfixture\tCDS\t13\t21\t.\t+\t0\t.");
		ProteinSupportCompare.main(new String[]{"mode=stops-v1", "fasta="+fasta, "ref="+ref, "baseline="+before, "candidate="+after, "out="+dir+"/paired"});
		final HashMap<String,String> metrics=metrics(dir+"/paired.summary.tsv");
		check("1".equals(metrics.get("reference.duplicate_coordinate_rows")), "Two identical rows contain one duplicate");
		check("1".equals(metrics.get("reference.primary_stops")), "Duplicate reference stop denominator is one");
		check("2".equals(metrics.get("candidate.primary_stops")) && "3".equals(metrics.get("candidate.primary_rows")), "Alternatives collapse only for stop denominators");
		check("0.5".equals(metrics.get("candidate.stop_reference_precision")) && "1.0".equals(metrics.get("candidate.stop_reference_recall")), "Literal distinct intersection oracle: 1/2 precision, 1/1 recall");
		check("1".equals(metrics.get("delta.gained_primary_stops")) && "0".equals(metrics.get("delta.lost_primary_stops"))
			&& "1".equals(metrics.get("delta.same_stop_changed_starts")), "Start-only replacements and true gained stops require different delta categories");
		check("2".equals(metrics.get("reference.table_rows.fallback_11")), "Fallback translation code use must be explicit even on duplicate rows");
		write(dir+"/empty.gff", "# no annotated CDS");
		ProteinSupportCompare.main(new String[]{"mode=stops-v1", "fasta="+fasta, "ref="+dir+"/empty.gff", "out="+dir+"/empty"});
		check("NA".equals(metrics(dir+"/empty.summary.tsv").get("reference.stop_reference_recall")), "Zero denominator must print NA");
		write(dir+"/bad.gff", "ctg\tfixture\tCDS\t1\t12");
		fails("truncated GFF", new Runnable(){public void run(){ProteinStopEvaluation.readGff(dir+"/bad.gff");}});
		for(int strand=0; strand<2; strand++){
			final String prefix=dir+"/native"+strand, sign=strand==0 ? "+" : "-";
			write(prefix+".fa", ">ctg\n"+(strand==0 ? "ATGAAANNN" : "NNNTTTCAT"));
			final String gff="ctg\tBBTools\tCDS\t1\t9\t.\t"+sign+"\t0\tCDS,fr0,start:ATG,stop:NNN";
			write(prefix+".before.gff", gff); write(prefix+".after.gff", gff+";bbgap_reason=synthetic_stop");
			final boolean original=GffLine.paseSource;
			check("BBTools".equals(ProteinStopEvaluation.readGff(prefix+".before.gff").get(0).source), "File reader must retain native producer provenance");
			check(GffLine.paseSource==original, "Opt-in reading must restore the legacy global source-parsing flag");
			ProteinSupportCompare.main(new String[]{"mode=stops-v1", "fasta="+prefix+".fa", "ref="+dir+"/empty.gff",
				"baseline="+prefix+".before.gff", "candidate="+prefix+".after.gff", "out="+prefix});
			final HashMap<String,String> synthetic=metrics(prefix+".summary.tsv");
			check("1".equals(synthetic.get("baseline.rows.synthetic_stop")) && "1".equals(synthetic.get("candidate.rows.synthetic_stop"))
				&& "0".equals(synthetic.get("candidate.primary_stops")), "File-roundtrip native/zero ambiguous boundaries must be excluded equally on both strands");
		}
		write(dir+"/circular.fa", ">ctg\nATGAAATAAATGAAATAAATGAAATAAATGAAATAA");
		write(dir+"/circular.gff", "ctg\tfixture\tregion\t1\t36\t.\t+\t.\tIs_circular=true\n"
			+"ctg\tfixture\tCDS\t28\t45\t.\t+\t0\tID=wrap;transl_table=4");
		ProteinSupportCompare.main(new String[]{"mode=stops-v1", "fasta="+dir+"/circular.fa", "ref="+dir+"/circular.gff", "out="+dir+"/circular"});
		final HashMap<String,String> circle=metrics(dir+"/circular.summary.tsv");
		check("1".equals(circle.get("reference.rows.circular_origin")) && "1".equals(circle.get("reference.table_rows.declared_4"))
			&& "2".equals(circle.get("reference.reference_mask_intervals")), "Excluded circular rows retain declarations and split masks in emitted inventory");
		containsLine(dir+"/circular.reference.tsv", "stops-v1\t2\tCDS\tctg\t0\t28\t45\t0\t4\t0\t0\tcircular_origin\tID=wrap;transl_table=4");
		containsLine(dir+"/circular.contigs.tsv", "ctg\t36\ttrue");
		System.err.println("STOP_EVALUATION_IO_PASS distinct_precision=0.5 distinct_recall=1.0 empty_denominator=NA");
	}
	private static void containsLine(String path, String expected){
		final ByteFile file=ByteFile.makeByteFile(path, true); boolean found=false;
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				found|=expected.equals(new String(line, StandardCharsets.UTF_8));
			}
		}finally{check(!file.close(), "Literal inventory read must finish");}
		check(found, "Missing literal inventory row in "+path+": "+expected);
	}
	private static HashMap<String,String> metrics(String path){
		final HashMap<String,String> result=new HashMap<String,String>();
		final ByteFile file=ByteFile.makeByteFile(path, true);
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				final String[] fields=new String(line, StandardCharsets.UTF_8).split("\t");
				check(fields.length==4, "Summary must retain four columns");
				result.put(fields[1]+"."+fields[2], fields[3]);
			}
		}finally{check(!file.close(), "Fixture summary read must finish");}
		return result;
	}
	private static void write(String path, String text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false); writer.start(); writer.println(text);
		check(!writer.poisonAndWait(), "Fixture must be written completely: "+path);
	}

	private static void one(String dna, String strand, int phase, String attrs, String expected){
		final ProteinStopEvaluation.Dataset data=normalize(dna, true, cds(1, dna.length(), strand, phase, attrs));
		check(data.rows.size()==1 && expected.equals(data.rows.get(0).reason), "Expected "+expected+" for "+dna+" "+strand+" "+attrs+", got "+data.rows.get(0).reason);
		check(data.primary.size()==(expected.equals("primary") ? 1 : 0), "Only primary rows may enter the denominator");
		if(expected.equals("primary") && attrs.contains("partial=" ) && !attrs.contains("00")){
			check(data.rows.get(0).five && !data.rows.get(0).three, "5-prime partial oracle must respect minus orientation");
		}
	}
	private static void duplicatesAndMasks(){
		final String dna="ATGATGAAATAA";
		final ProteinStopEvaluation.Dataset ref=normalize(dna, true,
			cds(1, 12, "+", 0, "ID=x"), cds(1, 12, "+", 0, "ID=x"), cds(4, 12, "+", 0, "ID=y"));
		check(ref.rawCds==3 && ref.primary.size()==1 && ref.all.size()==1, "Duplicate rows/alternative starts must count as ONE reference stop");
		final ProteinStopEvaluation.Dataset calls=normalize(dna, false, cds(1, 12, "+", 0, "."), cds(4, 12, "+", 0, "."));
		ProteinStopEvaluation.applyMasks(calls, ref);
		check(calls.primary.size()==1 && ref.primary.keySet().equals(calls.primary.keySet()), "Both precision and recall share a distinct intersection");
		final ProteinStopEvaluation.Dataset excluded=normalize(dna, true, cds(1, 12, "+", 0, "exception=unmodeled"));
		ProteinStopEvaluation.applyMasks(calls, excluded);
		check(calls.primary.isEmpty() && calls.rawCds==2 && calls.all.size()==1, "Masked calls remain in raw inventory but not primary denominator");
		final ProteinStopEvaluation.Dataset mixed=normalize(dna, true, cds(1, 12, "+", 0, "."), cds(4, 12, "+", 0, "pseudo=true"));
		final ProteinStopEvaluation.Dataset match=normalize(dna, false, cds(1, 12, "+", 0, "."));
		ProteinStopEvaluation.applyMasks(match, mixed);
		check(match.primary.size()==1, "Eligible exact stop match wins over an overlapping exception mask");
		final ProteinStopEvaluation.Dataset none=normalize(dna, true);
		check(none.primary.isEmpty() && none.rows.isEmpty(), "Empty reference has a zero denominator, not invented negative labels");
		final ProteinStopEvaluation.Dataset unannotated=normalize(dna, false, cds(1, 12, "+", 0, "."));
		ProteinStopEvaluation.applyMasks(unannotated, none);
		check(unannotated.primary.size()==1, "Unannotated calls stay in the precision denominator as unresolved");
	}
	private static void parentAndMultipart(){
		final String dna="ATGAAATAAATGAAATAA";
		final GffLine parent=line("gene", 1, 18, "+", -1, "ID=g;pseudo=true");
		final ProteinStopEvaluation.Dataset pseudo=normalize(dna, true, parent, cds(1, 9, "+", 0, "Parent=g;ID=x"));
		check(pseudo.primary.isEmpty() && pseudo.rows.get(0).reason.equals("pseudo") && pseudo.maskRows.size()==1, "Parent-only pseudo marks CDS and supplies a mask even without CDS");
		final ProteinStopEvaluation.Dataset parts=normalize(dna, true, cds(1, 9, "+", 0, "ID=x"), cds(10, 18, "+", 0, "ID=x"));
		check(parts.primary.isEmpty() && parts.rawCds==2 && parts.rows.get(0).reason.equals("multipart"), "Segment boundaries cannot count as ordinary biological stops");
		final ProteinStopEvaluation.Dataset noId=normalize(dna, true,
			line("gene", 1, 18, "+", -1, "ID=g"), cds(1, 9, "+", 0, "Parent=g"), cds(10, 18, "+", 0, "Parent=g"));
		check(noId.primary.isEmpty() && noId.rows.get(0).reason.equals("multipart"), "Missing CDS IDs must not hide multi-segment parent linkage");
		final ProteinStopEvaluation.Dataset loci=normalize(dna, true,
			line("gene", 1, 9, "+", -1, "ID=a"), line("gene", 10, 18, "+", -1, "ID=b"),
			cds(1, 9, "+", 0, "ID=x;Parent=a;protein_id=same"), cds(10, 18, "+", 0, "ID=x;Parent=b;protein_id=same"));
		check(loci.primary.size()==2, "Shared protein accession/CDS name across different parents cannot merge loci");
		// Literal geometry of the two archived pilot multipart proteins; no sequence assertion is made.
		final HashMap<String,byte[]> seq=new HashMap<String,byte[]>(); seq.put("ctg", new byte[408928]);
		final ArrayList<GffLine> pilot=new ArrayList<GffLine>();
		pilot.add(cds(377615, 377686, "-", 0, "ID=prfB;exception=ribosomal slippage"));
		pilot.add(cds(376591, 377613, "-", 0, "ID=prfB;exception=ribosomal slippage"));
		pilot.add(cds(407535, 407985, "+", 0, "ID=IS3;exception=ribosomal slippage"));
		pilot.add(cds(407985, 408928, "+", 2, "ID=IS3;exception=ribosomal slippage"));
		final ProteinStopEvaluation.Dataset archived=ProteinStopEvaluation.normalize(pilot, seq, 11, true);
		check(archived.rawCds==4 && archived.primary.isEmpty() && archived.masks.size()==2, "Four slippage segments are two excluded loci, not four primary reference stops");
	}
	private static void failures(){
		fails("bounds", new Runnable(){public void run(){normalize("ATGAAATAA", true, cds(1, 10, "+", 0, "."));}});
		fails("invalid predicted stop", new Runnable(){public void run(){normalize("ATGAAAAAA", false, cds(1, 9, "+", 0, "."));}});
		fails("invalid predicted frame", new Runnable(){public void run(){normalize("AATGAAATAA", false, cds(1, 10, "+", 0, "."));}});
		fails("partial syntax", new Runnable(){public void run(){normalize("ATGAAATAA", true, cds(1, 9, "+", 0, "partial=maybe"));}});
		fails("missing contig", new Runnable(){public void run(){
			final GffLine row=cds(1, 9, "+", 0, "."); row.seqid="missing"; normalize("ATGAAATAA", true, row);
		}});
	}
	private static void circularOrigin(){
		for(final String strand : new String[]{"+", "-"}){
			final String unit=strand.equals("+") ? "ATGAAATAA" : "TTATTTCAT", dna=unit+unit+unit+unit;
			final GffLine region=line("region", 1, 36, "+", -1, "Is_circular=true");
			final GffLine wrapped=cds(28, 45, strand, 0, "ID=w;transl_table=4");
			final ProteinStopEvaluation.Dataset ref=normalize(dna, true, region, wrapped);
			check(ref.rawCds==1 && ref.primary.isEmpty() && ref.masks.size()==2
				&& "circular_origin".equals(ref.rows.get(0).reason) && ref.rows.get(0).table==4,
				"Single-wrap reference interval must be excluded before sequence translation, preserving table");
			check(ref.rows.get(0).g.start==28 && ref.rows.get(0).g.stop==45
				&& ref.all.containsKey("ctg\t"+(strand.equals("+") ? "0\t45" : "1\t28")),
				"Raw coordinates and oriented virtual endpoint must remain unchanged on both strands");
			maskEndsOnly(dna, strand, ref);
			final ProteinStopEvaluation.Dataset override=normalize(dna, true, region, wrapped, cds(1, 9, strand, 0, "."));
			final ProteinStopEvaluation.Dataset matched=normalize(dna, false, cds(1, 9, strand, 0, "."));
			ProteinStopEvaluation.applyMasks(matched, override);
			check(matched.primary.size()==1, "Exact primary reference stop overrides a circular mask on either strand");
			final GffLine parent=line("gene", 28, 45, strand, -1, "ID=g");
			final ProteinStopEvaluation.Dataset children=normalize(dna, true, region, parent,
				cds(28, 36, strand, 0, "ID=c;Parent=g"), cds(1, 9, strand, 0, "ID=c;Parent=g"));
			check(children.rawCds==2 && children.masks.size()==2 && children.primary.isEmpty()
				&& "circular_parent".equals(children.rows.get(0).reason) && "circular_parent".equals(children.rows.get(1).reason),
				"Bounded children on both physical ends inherit explicit circular-parent exclusion, never a whole-contig hull");
			maskEndsOnly(dna, strand, children);
			final ProteinStopEvaluation.Dataset unstranded=normalize(dna, true, region,
				line("gene", 28, 45, ".", -1, "ID=u"), cds(1, 9, strand, 0, "Parent=u"));
			check("circular_parent".equals(unstranded.rows.get(0).reason), "Unknown parent strand remains compatible with known child strand");
			maskEndsOnly(dna, strand, unstranded);
			final ProteinStopEvaluation.Dataset chain=normalize(dna, true, region,
				line("mRNA", 28, 45, strand, -1, "ID=r;Parent=ancestor"), cds(1, 9, strand, 0, "Parent=r"));
			check("parent_chain".equals(chain.rows.get(0).reason) && chain.masks.size()==2, "Validated wrapping immediate parent retains parent-chain exclusion with bounded masks");
			maskEndsOnly(dna, strand, chain);
			final ProteinStopEvaluation.Dataset simple=normalize(dna, true, region, parent, cds(28, 45, strand, 0, "Parent=g"));
			check("circular_origin".equals(simple.rows.get(0).reason) && simple.masks.size()==2, "Matching virtual parent/CDS is the simple NCBI representation");
			final GffLine pseudo=line("gene", 28, 45, strand, -1, "ID=p;pseudo=true");
			final ProteinStopEvaluation.Dataset masks=normalize(dna, true, region, pseudo);
			check(masks.rawCds==0 && masks.maskRows.size()==1 && masks.masks.size()==2, "Wrapping pseudo non-CDS feature supplies two masks without inventing a CDS");
			maskEndsOnly(dna, strand, masks);
			final ProteinStopEvaluation.Dataset inherited=normalize(dna, true, region, pseudo, cds(28, 45, strand, 0, "Parent=p"));
			check("pseudo".equals(inherited.rows.get(0).reason), "Parent pseudo keeps precedence over circular-origin exclusion");
			final ProteinStopEvaluation.Dataset direct=normalize(dna, true, region, cds(28, 45, strand, 0, "pseudo=true"));
			check("pseudo".equals(direct.rows.get(0).reason), "CDS pseudo keeps precedence over circular-origin exclusion");
			final String invalid="Invalid circular-origin interval or missing validated region:";
			fails("predictions cannot wrap even with circular region", invalid, new Runnable(){public void run(){normalize(dna, false, region, wrapped);}});
			fails("missing circular marker", invalid, new Runnable(){public void run(){normalize(dna, true, wrapped);}});
			fails("multi-revolution interval", invalid, new Runnable(){public void run(){normalize(dna, true, region, cds(28, 65, strand, 0, "."));}});
			fails("large end must not overflow span", invalid, new Runnable(){public void run(){normalize(dna, true, region, cds(28, Integer.MAX_VALUE, strand, 0, "."));}});
			fails("virtual start outside actual sequence", invalid, new Runnable(){public void run(){normalize(dna, true, region, cds(37, 45, strand, 0, "."));}});
			fails("wrong marker feature", invalid, new Runnable(){public void run(){normalize(dna, true, line("gene", 1, 36, "+", -1, "Is_circular=true"), wrapped);}});
			fails("short region", "Circular region must match full FASTA length:", new Runnable(){public void run(){normalize(dna, true, line("region", 1, 35, "+", -1, "Is_circular=true"), wrapped);}});
			fails("long region", "Invalid 1-based inclusive GFF interval:", new Runnable(){public void run(){normalize(dna, true, line("region", 1, 37, "+", -1, "Is_circular=true"), wrapped);}});
			fails("offset region", "Circular region must match full FASTA length:", new Runnable(){public void run(){normalize(dna, true, line("region", 2, 36, "+", -1, "Is_circular=true"), wrapped);}});
			fails("conflicting regions", "Conflicting circular region markers:", new Runnable(){public void run(){normalize(dna, true, region, line("region", 1, 36, "+", -1, "Is_circular=false"), wrapped);}});
			fails("invalid parent before inconsistent-strand exclusion", invalid, new Runnable(){public void run(){normalize(dna, true, region,
				line("gene", 28, 65, strand.equals("+") ? "-" : "+", -1, "ID=bad"), cds(1, 9, strand, 0, "Parent=bad"));}});
		}
		System.err.println("STOP_EVALUATION_CIRCULAR_PASS strands=2 split_masks=true parent_wrap=true raw_coordinates=true strict_predictions=true");
	}
	private static void maskEndsOnly(String dna, String strand, ProteinStopEvaluation.Dataset ref){
		final ProteinStopEvaluation.Dataset calls=normalize(dna, false,
			cds(1, 9, strand, 0, "."), cds(10, 18, strand, 0, "."), cds(28, 36, strand, 0, "."));
		ProteinStopEvaluation.applyMasks(calls, ref);
		check("reference_exception_mask".equals(calls.rows.get(0).reason) && calls.rows.get(1).primary()
			&& "reference_exception_mask".equals(calls.rows.get(2).reason), "Both circular end pieces mask, but ordinary middle-of-contig stops must remain primary");
		// Opposite-strand terminal needs valid sequence, so use a literal excluded-free fixture row to test mask strand only.
		final GffLine opposite=cds(1, 9, strand.equals("+") ? "-" : "+", 0, ".");
		final ProteinStopEvaluation.Dataset other=new ProteinStopEvaluation.Dataset(false);
		final ProteinStopEvaluation.Row row=new ProteinStopEvaluation.Row(opposite, 1); row.reason="primary";
		other.rows.add(row); other.rawCds=1; other.index(); ProteinStopEvaluation.applyMasks(other, ref);
		check(other.primary.size()==1, "Circular masks must not leak onto the opposite strand");
	}
	private static void nativeSynthetic(){
		for(String strand : new String[]{"+", "-"}){
			final String dna=strand.equals("+") ? "ATGAAANNN" : "NNNTTTCAT";
			final GffLine before=cds(1, 9, strand, 0, "CDS,fr0,start:ATG,stop:NNN"); before.source="BBTools";
			final GffLine after=cds(1, 9, strand, 0, "CDS,fr0,start:ATG,stop:NNN;bbgap_reason=synthetic_stop"); after.source="BBTools";
			final ProteinStopEvaluation.Dataset a=normalize(dna, false, before), b=normalize(dna, false, after);
			check(a.primary.isEmpty() && b.primary.isEmpty() && a.rows.get(0).reason.equals("synthetic_stop")
				&& b.rows.get(0).reason.equals("synthetic_stop"), "Native and annotated ambiguous stops must classify equally on both strands");
		}
		fails("false unknown-stop label", new Runnable(){public void run(){
			final GffLine g=cds(1, 9, "+", 0, "CDS,fr0,start:ATG,stop:NNN"); g.source="BBTools";
			normalize("ATGAAATAA", false, g);
		}});
		fails("unknown-stop string in a note", new Runnable(){public void run(){
			final GffLine g=cds(1, 9, "+", 0, "Note=stop:NNN"); g.source="BBTools";
			normalize("ATGAAANNN", false, g);
		}});
		fails("unknown-stop label from unrelated producer", new Runnable(){public void run(){
			normalize("ATGAAANNN", false, cds(1, 9, "+", 0, "CDS,fr0,start:ATG,stop:NNN"));
		}});
	}
	private static ProteinStopEvaluation.Dataset normalize(String dna, boolean reference, GffLine... rows){
		final ArrayList<GffLine> list=new ArrayList<GffLine>(); for(GffLine row : rows){list.add(row);}
		final HashMap<String,byte[]> seq=new HashMap<String,byte[]>(); seq.put("ctg", dna.getBytes(StandardCharsets.US_ASCII));
		return ProteinStopEvaluation.normalize(list, seq, 11, reference);
	}
	private static GffLine cds(int start, int stop, String strand, int phase, String attrs){return line("CDS", start, stop, strand, phase, attrs);}
	private static GffLine line(String type, int start, int stop, String strand, int phase, String attrs){
		return new GffLine(("ctg\tfixture\t"+type+"\t"+start+"\t"+stop+"\t.\t"+strand+"\t"+(phase<0 ? "." : phase)+"\t"+attrs).getBytes(StandardCharsets.UTF_8));
	}
	private static void fails(String label, Runnable action){
		boolean failed=false; try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Must fail loudly: "+label);
	}
	private static void fails(String label, String diagnostic, Runnable action){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(diagnostic), "Wrong diagnostic for "+label+": "+expected.getMessage());
			failed=true;
		}
		check(failed, "Must fail loudly: "+label);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private ProteinStopEvaluationTest(){}
}
