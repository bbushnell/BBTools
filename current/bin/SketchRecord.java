package bin;

import json.JsonObject;
import shared.LineParserS2;
import structures.ByteBuilder;

/**
 * Represents a genomic sketch record with metadata for similarity analysis.
 * Stores sketch data including ANI, completeness, contamination, and taxonomic information.
 * Extends BinObject to inherit k-mer processing and taxonomic utilities.
 *
 * @author Brian Bushnell
 * @date Feb 4, 2025
 */
public class SketchRecord extends BinObject {
	
	/**
	 * Constructs a SketchRecord from a JsonObject containing sketch hit data.
	 * Extracts and sets all sketch metadata fields from the JSON object.
	 * @param hit JsonObject containing sketch hit data with keys: Matches, ANI, Complt, Contam, TaxID, taxName
	 */
	public SketchRecord(JsonObject hit) {
		setFrom(hit);
	}
	
	/**
	 * Populates this sketch record's fields from a JsonObject.
	 * Extracts numeric values for matches, ANI, completeness, contamination, and taxonomy ID.
	 * Sets genus taxonomy ID to -1 as placeholder for future tree-based lookup.
	 * @param hit JsonObject containing sketch metadata
	 */
	public void setFrom(JsonObject hit) {
		matches=hit.getLong("Matches").intValue();
		ani=hit.getDouble("ANI").floatValue();
		completeness=hit.getDouble("Complt").floatValue();
		contam=hit.getDouble("Contam").floatValue();
		taxid=hit.getLong("TaxID").intValue();
		taxName=hit.getString("taxName");
		genusTaxid=-1;//TODO, get from tree
	}
	
	/**
	 * Converts this sketch record to a ByteBuilder representation.
	 * Creates a new ByteBuilder and appends all sketch metadata to it.
	 * @return ByteBuilder containing formatted sketch record data
	 */
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		return appendTo(bb);
	}
	
	/**
	 * Appends this sketch record's data to an existing ByteBuilder.
	 * Formats all metadata fields with tab separation and 2 decimal places for float values.
	 * Applies name shrinking to taxonomic names longer than 40 characters.
	 *
	 * @param bb ByteBuilder to append data to
	 * @return The same ByteBuilder with appended sketch data
	 */
	public ByteBuilder appendTo(ByteBuilder bb) {
		bb.append("ANI: ").append(ani, 2);
		bb.tab().append("Complt: ").append(completeness, 2);
		bb.tab().append("Contam: ").append(contam, 2);
		bb.tab().append("Matches: ").append(matches);
		bb.tab().append("TaxID: ").append(taxid);
		bb.tab().append("Name: ").append(shrink(taxName));
		return bb;
	}
	
	/**
	 * Shortens taxonomic names longer than 40 characters to improve readability.
	 * Takes first 4 space-separated words and appends "..." if more words exist.
	 * Returns original name if shortened version is not actually shorter.
	 *
	 * @param n Taxonomic name to potentially shorten
	 * @return Shortened name or original if shortening provides no benefit
	 */
	private static String shrink(String n) {
		if(n==null || n.length()<40) {return n;}
		ByteBuilder bb=new ByteBuilder();
		LineParserS2 lp=new LineParserS2(' ');
		lp.set(n);
		for(int i=0; i<4 && lp.hasMore(); i++) {
			if(bb.length>0) {bb.space();}
			bb.append(lp.parseString());
		}
		if(lp.hasMore()) {bb.append("...");}
		return (bb.length()<n.length() ? bb.toString() : n);
	}

	/** Number of k-mer matches between query and reference sketch */
	int matches=-1;
	/** Completeness score indicating genome assembly completeness percentage */
	float completeness=-1;
	/** Contamination score indicating percentage of contaminating sequences */
	float contam=-1;
	/** Average Nucleotide Identity percentage between query and reference */
	float ani=-1;
	/** Taxonomic identifier for the sketch record */
	int taxid=-1;
	/** Genus-level taxonomic identifier, currently unused (-1) */
	int genusTaxid=-1;
	/** Human-readable taxonomic name associated with this sketch record */
	String taxName=null;
	
}
