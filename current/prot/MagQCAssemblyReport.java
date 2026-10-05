package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import bin.BinStats;
import stream.Read;
import structures.ByteBuilder;

/** Assembly statistics for public reporting, independent of the model input layout.
 * @author Yoimiya */
final class MagQCAssemblyReport {

	/** Copies existing sufficient statistics and computes Nx/Lx on whole FASTA records. */
	MagQCAssemblyReport(ArrayList<Read> reads, MagQCPreparedBin bin){
		this(reads, bin, false);
	}

	/** BBTools uses N for count and L for length; swapnl changes labels as in AssemblyStats2. */
	MagQCAssemblyReport(ArrayList<Read> reads, MagQCPreparedBin bin, boolean swapNL_){
		swapNL=swapNL_;
		assert(reads!=null && !reads.isEmpty() && bin!=null) :
			"Reporting follows successful assembly loading and native prepared-bin construction";
		contigs=reads.size(); length=bin.length;
		final int[] lengths=new int[contigs];
		long total=0;
		for(int i=0; i<contigs; i++){
			lengths[i]=reads.get(i).length(); total+=lengths[i];
		}
		if(total!=length){throw new IllegalStateException("Report and prepared-bin assembly lengths differ");}
		Arrays.sort(lengths);
		final long half=(length+1)/2, ninety=length-length/10;
		long sum=0;
		int count=0;
		for(int i=lengths.length-1; i>=0; i--){
			sum+=lengths[i]; count++;
			if(n50==0 && sum>=half){n50=count; l50=lengths[i];}
			if(sum>=ninety){n90=count; l90=lengths[i]; break;}
		}
		assert(n50>0 && n90>=n50) : "Positive assembly length must reach both cumulative base thresholds in descending contig order";
		final MagQCVectorMaker.Agg a=bin.agg;
		gc=a.acgt==0 ? Double.NaN : a.gc/(double)a.acgt;
		coding=a.coding/(double)length;
		cds=a.cds; r16=a.r16; r23=a.r23; r5=a.r5; trna=a.trna;
	}

	/** Appends fraction-valued statistics; missing classifier metadata is explicit NA. */
	void append(ByteBuilder out, MagQCAssemblyInput.Taxonomy tax, float comp, float contam){
		assert(out!=null && tax!=null) : "The completed bin owns its report buffer and classifier result";
		out.tab().append(tax.name).tab();
		if(tax.taxId<0){out.append("NA");}else{out.append(tax.taxId);}
		out.tab();
		if(tax.ani<0){out.append("NA");}else{out.appendSlow(tax.ani);}
		out.tab().append(contigs).tab().append(length).tab().append(n50).tab().append(l50);
		out.tab().append(n90).tab().append(l90).tab();
		if(Double.isNaN(gc)){out.append("NA");}else{out.appendSlow(gc);}
		out.tab().append(cds).tab().appendSlow(coding).tab().append(r16).tab().append(r23);
		out.tab().append(r5).tab().append(trna).tab().append(quality(comp, contam)).tab().append("NA");
	}

	/** Uses GradeBins' existing RNA-aware extended MIMAG classifier without clipping predictions. */
	String quality(float comp, float contam){
		if(!Float.isFinite(comp) || !Float.isFinite(contam)){
			throw new IllegalArgumentException("Quality classification requires finite predictions");
		}
		return BinStats.type(comp, contam, r16, r23, r5, trna, true);
	}

	/** Formats one assembly after success; errors are predicted absolute errors, not confidence intervals. */
	String human(MagQCAssemblyInput.Taxonomy tax, float comp, float contam, double compError, double contamError){
		assert(tax!=null && compError>=0 && contamError>=0) :
			"The public scorer validates taxonomy and scaled error heads before formatting";
		final StringBuilder out=new StringBuilder(512);
		line(out, "Name:", tax.name); line(out, "Phylum:", tax.phylum);
		line(out, "TaxID:", tax.taxId<0 ? "NA" : Long.toString(tax.taxId));
		line(out, "Size:", length+" bp"); line(out, "Contigs:", Integer.toString(contigs));
		line(out, "ANI:", tax.ani<0 ? "NA" : decimal(100*tax.ani, 2));
		line(out, swapNL ? "L50/N50:" : "N50/L50:", n50+"/"+l50);
		line(out, swapNL ? "L90/N90:" : "N90/L90:", n90+"/"+l90);
		line(out, "Completeness:", decimal(100.0*comp, 2)+" +-"+decimal(100*compError, 2));
		line(out, "Contamination:", decimal(100.0*contam, 2)+" +-"+decimal(100*contamError, 2));
		line(out, "GC:", Double.isNaN(gc) ? "NA" : decimal(gc, 3));
		line(out, "CDS:", Integer.toString(cds)); line(out, "Coding Density:", decimal(100*coding, 2));
		line(out, "16S:", Integer.toString(r16)); line(out, "23S:", Integer.toString(r23));
		line(out, "5S:", Integer.toString(r5)); line(out, "tRNA:", Integer.toString(trna));
		line(out, "18S:", "NA");//The serving caller disables 18S; zero would falsely claim absence.
		line(out, "Quality:", quality(comp, contam));
		return out.toString();
	}

	/** Aligns labels in a single small human-facing block. */
	private static void line(StringBuilder out, String label, String value){
		assert(label.length()<17) : "Report labels must fit the fixed human-readable label column";
		out.append(label);
		for(int i=label.length(); i<17; i++){out.append(' ');}
		out.append(value).append('\n');
	}

	/** Uses a stable decimal point regardless of the host locale. */
	private static String decimal(double value, int places){
		assert(Double.isFinite(value)) : "Missing metrics must be rendered as NA before decimal formatting";
		return String.format(Locale.ROOT, "%."+places+"f", value);
	}

	/** Keeps count before length in both conventions; unit-bearing names identify each TSV value. */
	static String columns(boolean swapNL){
		return "\treference_name\treference_taxid\tani_fraction\tcontigs\tgenome_size_bp"
			+(swapNL ? "\tl50_contigs\tn50_bp\tl90_contigs\tn90_bp" : "\tn50_contigs\tl50_bp\tn90_contigs\tl90_bp")
			+"\tgc_fraction\tcds\tcoding_density_fraction\tr16\tr23\tr5\ttrna\tquality\tr18";
	}
	private final boolean swapNL;
	private final int contigs, cds, r16, r23, r5, trna;
	private final long length;
	private final double gc, coding;
	int n50, l50, n90, l90;
}
