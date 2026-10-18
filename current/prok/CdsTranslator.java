package prok;

import java.util.Arrays;

import dna.AminoAcid;
import dna.GeneticCode;
import gff.GffLine;

/** Translation of contiguous GFF CDS intervals without changing their annotation.
 * Binary partial flags describe low/high genomic edges, not biological 5'/3'.
 * Unsupported biological forms return a reason; invalid coordinates fail loudly.
 * @author Keqing
 */
public final class CdsTranslator {

	public static Translation translate(final GffLine feature, final byte[] scaffold, final int fallbackTable){
		if(feature==null || scaffold==null || feature.start<1 || feature.stop<feature.start || feature.stop>scaffold.length){
			throw new IllegalArgumentException("CDS coordinates must be 1-based inclusive and inside the supplied scaffold");
		}
		final String declared=attribute(feature.attributes, "transl_table");
		final int table=declared==null ? fallbackTable : Integer.parseInt(declared);
		final String partial=attribute(feature.attributes, "partial");
		//Use the existing GFF reader's duplicate/value validation before interpreting sides.
		feature.partial();
		final boolean left=partial!=null && (partial.equals("true") || partial.charAt(0)=='1');
		final boolean right=partial!=null && (partial.equals("true") || (partial.length()==2 && partial.charAt(1)=='1'));
		final boolean five=feature.strand==GffLine.MINUS ? right : left;
		final boolean three=feature.strand==GffLine.MINUS ? left : right;
		if(table!=4 && table!=11 && table!=25){return unsupported(table, five, three, "unsupported_translation_table");}
		if(feature.strand!=GffLine.PLUS && feature.strand!=GffLine.MINUS){return unsupported(table, five, three, "unknown_strand");}
		if(feature.phase<0 || feature.phase>2){return unsupported(table, five, three, "unknown_phase");}
		if(attribute(feature.attributes, "transl_except")!=null || attribute(feature.attributes, "exception")!=null){
			return unsupported(table, five, three, "translation_exception");
		}
		final byte[] bases=Arrays.copyOfRange(scaffold, feature.start-1, feature.stop);
		if(feature.strand==GffLine.MINUS){AminoAcid.reverseComplementBasesInPlace(bases);}
		final int length=bases.length-feature.phase;
		if(length<3){return unsupported(table, five, three, "empty_translation");}
		if(length%3!=0 && !three){return unsupported(table, five, three, "incomplete_terminal_codon");}
		final GeneticCode code=GeneticCode.forTable(table);
		int count=length/3;
		if(!three){
			if(!code.isStop(GeneticCode.codon(bases, feature.phase+3*(count-1)))){
				return unsupported(table, five, three, "missing_terminal_stop");
			}
			count--;
		}
		if(count<1){return unsupported(table, five, three, "empty_translation");}
		final byte[] peptide=new byte[count];
		for(int i=0; i<count; i++){
			final int codon=GeneticCode.codon(bases, feature.phase+3*i);
			if(codon<0){return unsupported(table, five, three, "ambiguous_codon");}
			if(code.isStop(codon)){
				return unsupported(table, five, three, i==count-1 && three ? "stop_at_partial_end" : "internal_stop");
			}
			final boolean initiation=i==0 && !five;
			if(initiation && !code.isStart(codon)){return unsupported(table, five, three, "invalid_start_codon");}
			peptide[i]=code.translateCodon(codon, initiation);
		}
		assert(peptide.length>0) : "Protein-family alignment requires a nonempty peptide after terminal-stop removal";
		return new Translation(peptide, table, five, three, null);
	}

	/** Exact semicolon-delimited key lookup; duplicate recognized keys are invalid. */
	static String attribute(final String attributes, final String key){
		if(attributes==null || attributes.equals(".")){return null;}
		String value=null;
		for(int start=0; start<attributes.length(); ){
			int end=attributes.indexOf(';', start);
			if(end<0){end=attributes.length();}
			final int eq=attributes.indexOf('=', start);
			if(eq>=start && eq<end && eq-start==key.length() && attributes.regionMatches(start, key, 0, key.length())){
				if(value!=null){throw new IllegalArgumentException("Duplicate GFF attribute: "+key);}
				value=attributes.substring(eq+1, end);
				if(value.isEmpty()){throw new IllegalArgumentException("Empty GFF attribute: "+key);}
			}
			start=end+1;
		}
		return value;
	}

	private static Translation unsupported(int table, boolean five, boolean three, String reason){
		assert(reason!=null && !reason.isEmpty()) : "Unsupported CDSs must retain a diagnostic reason in annotation output";
		return new Translation(null, table, five, three, reason);
	}

	public static final class Translation {
		Translation(byte[] peptide_, int table_, boolean five_, boolean three_, String reason_){
			peptide=peptide_; table=table_; partialFive=five_; partialThree=three_; reason=reason_;
		}
		public final byte[] peptide;
		final int table;
		final boolean partialFive, partialThree;
		public final String reason;
	}
	private CdsTranslator(){}
}
