package prot;

import java.util.ArrayList;
import java.util.Arrays;

import parse.LineParser1;
import stream.Read;
import structures.ByteBuilder;

/** Reporting fixtures for uneven contigs, unavailable metadata, fractions and RNA quality gates.
 * @author Yoimiya */
public final class MagQCAssemblyReportTest {

	/** Exercises the same formatter used by single-bin stderr and every assembly TSV row. */
	public static void main(String[] args){
		final ArrayList<Read> reads=new ArrayList<Read>();
		for(int length:new int[]{1, 5, 4}){
			final byte[] bases=new byte[length]; Arrays.fill(bases, (byte)'A');
			reads.add(new Read(bases, null, "contig"+length, length));
		}
		final MagQCPreparedBin bin=new MagQCPreparedBin(1);
		bin.length=10; bin.agg.gc=3; bin.agg.acgt=8; bin.agg.coding=11; bin.agg.cds=4;
		bin.agg.r16=2; bin.agg.r23=1; bin.agg.r5=3; bin.agg.trna=18;
		final MagQCAssemblyReport report=new MagQCAssemblyReport(reads, bin);
		check(report.n50==5 && report.l50==1 && report.n90==4 && report.l90==2,
			"Nx uses descending whole-contig lengths, Lx counts, and inclusive cumulative boundaries");
		final MagQCAssemblyInput.Taxonomy tax=new MagQCAssemblyInput.Taxonomy("Bacteria", "Bacillota", "reference", 123, .8351);
		final ByteBuilder row=new ByteBuilder(); report.append(row, tax, .997f, .006f);
		final LineParser1 p=new LineParser1('\t'); p.set(row.toBytes());
		check(p.terms()==18 && p.termEquals("0.8351", 3) && p.termEquals("0.375", 10) && p.termEquals("1.1", 12),
			"TSV must retain fraction units, use GC/ACGT, and preserve overlapping coding bases: "+row);
		check(p.termEquals("UHQ", 17), "RNA-complete high-quality fixture must use existing extended tier");
		final String human=report.human(tax, .997f, .006f, .0051, .0033);
		check(human.contains("83.51") && human.contains("99.70 +-0.51") && human.contains("0.60 +-0.33") && human.contains("110.00"),
			"Human percentages must scale values and error heads exactly once");
		bin.agg.trna=17; bin.agg.acgt=0; bin.agg.gc=0;
		final MagQCAssemblyReport missing=new MagQCAssemblyReport(reads, bin);
		check(missing.quality(.997f, .006f).equals("MQ"), "Missing RNA support must prevent HQ assignment");
		row.clear(); missing.append(row, new MagQCAssemblyInput.Taxonomy("unknown", "unknown"), .997f, .006f);
		p.set(row.toBytes());
		check(p.termEquals("NA", 1) && p.termEquals("NA", 2) && p.termEquals("NA", 3) && p.termEquals("NA", 10),
			"Unknown reference and undefined GC must not be invented as zero");
		final String prefix="#Query1\nmagqc_bin\t0.5\t10\t3\treference\t123\t0.5\t20\t4\tspecies"
			+"\t0\t0\t0\t0\t0\t0\t0";
		final String suffix="\td__Bacteria;p__Bacillota\tdomain\t.\n";
		for(String ssu:new String[]{"", "\t0.7"}){
			final MagQCAssemblyInput.Taxonomy parsed=MagQCAssemblyInput.parseResponse(
				prefix+ssu+"\t0.8351\t0.5\t0.4\t0.3\t10\t20"+suffix, 10, 3);
			check(parsed.name.equals("reference") && parsed.taxId==123 && parsed.ani==.8351,
				"ANI must be the first of six DDL metrics regardless of the optional SSU column");
			check(MagQCAssemblyInput.parseResponse(prefix+ssu+suffix, 10, 3).ani<0,
				"A response without DDL must report unavailable ANI");
		}
		System.out.println("MagQCAssemblyReportTest PASS: Nx/Lx, fractions, taxonomy metadata, missing values and RNA quality");
	}

	/** Keeps output-contract checks active even when assertions are disabled for malformed-input tests. */
	private static void check(boolean pass, String message){if(!pass){throw new AssertionError(message);}}
}
