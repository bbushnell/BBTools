package prot;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import shared.Timer;
import structures.IntList;

/**
 * Command-line front end for {@link MagQC}: reads a marker-vector TSV (as written
 * by {@code markervector.sh} / {@link MarkerVectorCLI}), reconstructs the
 * {@link MarkerVector}, runs the CheckM1-style completeness/contamination oracle,
 * and writes the full QC report with every frozen field (§4d).
 *
 * <p>This is a thin CLI/testing wrapper; the estimation lives in {@link MagQC} so
 * it can be called in-process from other BBTools code with no disk round-trip.
 * The shipped estimator is the §3d neural net — this reports the oracle it is
 * validated against.</p>
 *
 * <p>Usage: {@code magqc.sh vector=<vec.tsv> out=<report.tsv>}</p>
 *
 * @author Eru
 */
public final class MagQCCLI {

	/**
	 * Program entry point: parses arguments, reads the vector TSV, runs the oracle,
	 * writes the report.
	 * @param args Command-line arguments (flag=value).
	 */
	public static void main(String[] args){
		final Timer t=new Timer();
		if(args.length==0){printUsageAndExit();}

		String vectorFile=null, out=null;
		boolean overwrite=true;

		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("vector") || a.equals("vec") || a.equals("in") || a.equals("i")){vectorFile=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(vectorFile==null){
			System.err.println("Error: vector= is required.\n");
			printUsageAndExit();
		}

		final MarkerVector mv=readVector(vectorFile);
		System.err.println("Loaded vector dim="+mv.dimension()+" from "+vectorFile);

		final MagQCResult r=new MagQC().estimate(mv, null);
		writeReport(r, out, overwrite);

		t.stop();
		System.err.println("Completeness="+fmt(r.completeness)+"% contamination="+
			fmt(r.contamination)+"% (excess-copy) in "+t);
	}

	/**
	 * Reads a marker-vector TSV and reconstructs the {@link MarkerVector}. Data rows
	 * are {@code family_id \t representative \t count}; {@code #}-prefixed lines are
	 * the header and derived scalars (domain / proteins matched-unmatched are read
	 * back, the rest are recomputed).
	 * @param vectorFile Path to the vector TSV.
	 * @return The reconstructed vector.
	 */
	static MarkerVector readVector(final String vectorFile){
		final IntList familyIds=new IntList();
		final ArrayList<String> reps=new ArrayList<String>();
		final IntList counts=new IntList();
		String domain="NA";
		int matched=-1, unmatched=-1;

		final ByteFile bf=ByteFile.makeByteFile(vectorFile, false);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				lp.set(line);
				if(lp.terms()>=2){
					final String key=lp.parseString(0);
					if(key.equals("#domain")){domain=lp.parseString(1);}
					else if(key.equals("#proteins_matched")){matched=lp.parseInt(1);}
					else if(key.equals("#proteins_unmatched")){unmatched=lp.parseInt(1);}
				}
				continue;//header + all scalars
			}
			lp.set(line);
			if(lp.terms()<3){
				throw new RuntimeException("Malformed vector row (need 3 columns): '"+new String(line)+"'");
			}
			familyIds.add(lp.parseInt(0));
			reps.add(lp.parseString(1));
			counts.add(lp.parseInt(2));
		}
		bf.close();

		if(counts.size==0){
			throw new RuntimeException("No vector rows in "+vectorFile+
				" (expected 'family_id<tab>representative<tab>count' rows).");
		}

		final int n=counts.size;
		final int[] c=new int[n];
		final int[] fid=new int[n];
		final String[] rep=new String[n];
		for(int i=0; i<n; i++){
			c[i]=counts.get(i);
			fid[i]=familyIds.get(i);
			rep[i]=reps.get(i);
		}
		//matched/unmatched are context only; default to 0 when the scalars were absent.
		return new MarkerVector(c, fid, rep, domain,
			(matched<0 ? 0 : matched), (unmatched<0 ? 0 : unmatched));
	}

	/** Writes the QC report to stdout or a file. */
	static void writeReport(final MagQCResult r, final String out, final boolean overwrite){
		final String[] lines=reportLines(r);
		if(out==null || out.equalsIgnoreCase("stdout")){
			for(final String s : lines){System.out.println(s);}
			return;
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TEXT, null, false,
			overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(final String s : lines){bsw.println(s);}
		bsw.poisonAndWait();
	}

	/**
	 * Builds the report as {@code field \t value} lines — the two percentages plus
	 * every field the freeze requires to travel with them (§4d).
	 * @param r The QC result.
	 * @return Report lines.
	 */
	static String[] reportLines(final MagQCResult r){
		return new String[]{
			"#MAG-QC completeness/contamination (CheckM1-style oracle)",
			"metric\tvalue",
			"completeness_pct\t"+fmt(r.completeness),
			"contamination_pct\t"+fmt(r.contamination),
			"contamination_pct_multicopy\t"+fmt(r.contaminationMultiCopy),
			"expected_markers\t"+r.expectedMarkers,
			"detected_markers\t"+r.detectedMarkers,
			"multicopy_markers\t"+r.multiCopyMarkers,
			"excess_copies\t"+r.excessCopies,
			"effective_denominator\t"+r.effectiveDenominator,
			"domain_assignment\t"+r.domainAssignment,
			"marker_set_id\t"+r.markerSetId,
			"lineage_taxid\t"+r.lineageTaxID,
			"rank\t"+r.rank,
			"assignment_confidence\t"+fmt(r.assignmentConfidence),
			"assignment_confidence_model\t"+r.assignmentConfidenceModel,
			"ood_status\t"+r.oodStatus,
			"sufficient_evidence\t"+r.sufficientEvidence
		};
	}

	/** Formats a percentage/confidence to 4 decimals, or "NA" for NaN. */
	static String fmt(final double d){
		if(Double.isNaN(d)){return "NA";}
		return String.format("%.4f", d);
	}

	/** Prints usage text and exits. */
	static void printUsageAndExit(){
		System.err.println(
			"MagQC (BBTools prot package) — marker vector -> completeness/contamination\n"+
			"\n"+
			"Computes the CheckM1-style marker-counting completeness/contamination estimate\n"+
			"(the ORACLE the MAG-QC neural net is validated against) from a marker-vector\n"+
			"TSV written by markervector.sh.\n"+
			"\n"+
			"Usage: magqc.sh vector=<vec.tsv> out=<report.tsv>\n"+
			"\n"+
			"Required:\n"+
			"  vector=   Marker-vector TSV from markervector.sh (family_id/rep/count rows).\n"+
			"Optional:\n"+
			"  out=      Output report TSV. Default: stdout.\n"+
			"  ow=       Overwrite output (t/f, default t).\n"+
			"\n"+
			"Output: 'metric<tab>value' lines — completeness_pct, contamination_pct (headline,\n"+
			"excess-copy), contamination_pct_multicopy (secondary), the raw counts\n"+
			"(expected/detected/multicopy markers, excess copies, effective denominator),\n"+
			"domain_assignment, marker_set_id, lineage_taxid, rank, assignment_confidence,\n"+
			"ood_status (fixed 'unknown' until OD-9), and sufficient_evidence.\n"+
			"Note: the shipped estimator is a neural net; this is the counting oracle.\n");
		System.exit(0);
	}
}
