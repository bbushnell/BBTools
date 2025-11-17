package icecream;

import stream.Read;
import structures.ByteBuilder;

/**
 * Constructs and manipulates PacBio sequencing read metadata with multi-subread support.
 * Parses and reconstructs complex sequencing read information, handling movie start/stop
 * coordinates, number of passes, subreads, and error rates. Supports combining multiple
 * read builders and generating formatted read headers.
 *
 * @author Brian Bushnell
 */
class ReadBuilder {
	
	/**
	 * Constructs a ReadBuilder with byte array bases and metadata.
	 * @param bases_ Sequence bases as byte array
	 * @param passes_ Number of sequencing passes
	 * @param movieStart_ Starting position in movie coordinates
	 * @param zmw_ Zero-Mode Waveguide identifier
	 */
	public ReadBuilder(byte[] bases_, float passes_, int movieStart_, long zmw_) {
		this(new ByteBuilder(bases_), passes_, movieStart_, zmw_);
	}
	
	/**
	 * Constructs a ReadBuilder with ByteBuilder and metadata.
	 * Sets movie stop position based on sequence length.
	 * @param bases_ Sequence bases as ByteBuilder
	 * @param passes_ Number of sequencing passes
	 * @param movieStart_ Starting position in movie coordinates
	 * @param zmw_ Zero-Mode Waveguide identifier
	 */
	public ReadBuilder(ByteBuilder bases_, float passes_, int movieStart_, long zmw_) {
		bases=bases_;
		passes=passes_;
		movieStart=movieStart_;
		movieStop=movieStart+bases.length();
		zmw=zmw_;

		fullPasses=passes<1 ? 0 : 1;
	}
	
	/**
	 * Determines if a read ID represents an "ice cream" read with multiple subreads.
	 * Parses tab-delimited ID to extract subread count.
	 * @param id Read identifier string with tab-delimited metadata
	 * @return true if subreads > 1, indicating ice cream contamination
	 */
	public static boolean isIceCream(String id){
		String[] terms=id.split("\t");
		int subreads=Integer.parseInt(terms[3].split("=")[1]);
		return subreads>1;
	}
	
	/**
	 * Parses a Read object to extract PacBio metadata and create ReadBuilder.
	 * Extracts movie coordinates, ZMW ID, passes, subreads, adapters, and error rate
	 * from tab-delimited read ID format.
	 * @param r Read object with formatted PacBio ID
	 * @return ReadBuilder populated with parsed metadata
	 */
	public static ReadBuilder parse(Read r) {
		ByteBuilder bases=new ByteBuilder(r.bases);
		String[] terms=r.id.split("\t");
		String[] name=terms[0].split("/");
		String[] position=name[2].split("_");
		
		int movieStart=Integer.parseInt(position[0]);
		int movieStop=Integer.parseInt(position[1]);
		long zmw=Long.parseLong(name[1]);
		
		float passes=Float.parseFloat(terms[1].split("=")[1]);
		int fullPasses=Integer.parseInt(terms[2].split("=")[1]);
		int subreads=Integer.parseInt(terms[3].split("=")[1]);
		int missing=Integer.parseInt(terms[4].split("=")[1]);
		int adapters=Integer.parseInt(terms[5].split("=")[1]);
		float errorRate=(terms.length<7 ? 0 : Float.parseFloat(terms[6].split("=")[1]));
		
		ReadBuilder rb=new ReadBuilder(bases, passes, movieStart, zmw);
		rb.movieStop=movieStop;
		rb.passes=passes;
		rb.fullPasses=fullPasses;
		rb.subreads=subreads;
		rb.missing=missing;
		rb.adapters=adapters;
		rb.errorRate=errorRate;
		return rb;
	}
	
	@Override
	public String toString(){
		return toHeader().toString();
	}
	
	/**
	 * Generates formatted PacBio read header with all metadata fields.
	 * Format: m1_2_3/zmw/start_stop passes=X fullPasses=X subreads=X missing=X
	 * adapters=X errorRate=X
	 * @return ByteBuilder containing formatted header
	 */
	public ByteBuilder toHeader(){
		ByteBuilder id=new ByteBuilder(200);
		id.append("m1_2_3/");
		id.append(zmw).append('/').append(movieStart).append('_').append(movieStop);
		id.tab().append("passes=").append(passes, 2);
		id.tab().append("fullPasses=").append(fullPasses);
		id.tab().append("subreads=").append(subreads);
		id.tab().append("missing=").append(missing);
		id.tab().append("adapters=").append(adapters);
		id.tab().append("errorRate=").append(errorRate, 3);
		return id;
	}
	
	/** Returns the length of the sequence bases.
	 * @return Number of bases in the sequence */
	public int length() {
		return bases.length();
	}
	
	/**
	 * Combines this ReadBuilder with another by appending bases and summing metadata.
	 * Updates movie stop position, missing bases, adapters, full passes, subreads,
	 * and total passes.
	 * @param rb ReadBuilder to append to this one
	 */
	void add(ReadBuilder rb){
		bases.append(rb.bases);
		
		movieStop+=rb.length();
		missing+=rb.missing;
		adapters+=rb.adapters;
		fullPasses+=rb.fullPasses;
		subreads+=rb.subreads;
		passes+=rb.passes;
	}
	
	/**
	 * Converts ReadBuilder to a Read object with formatted header.
	 * Creates Read with bases, null quality scores, and complete metadata header.
	 * @return Read object ready for downstream processing
	 */
	Read toRead() {
		//Example: m54283_190403_183820/4194374/919_2614
		//Run ID is m54283_190403_183820
		//zmw ID is 4194374.
		//Read start/stop coordinates are 919_2614
		
		ByteBuilder id=toHeader();
		Read r=new Read(bases.toBytes(), null, id.toString(), 0);
		return r;
	}
	
	/** Sequence bases stored as ByteBuilder for efficient manipulation */
	ByteBuilder bases;

	/** Zero-Mode Waveguide identifier for this read */
	final long zmw;
	/** Starting coordinate position in the movie */
	final int movieStart;
	/** Ending coordinate position in the movie */
	int movieStop;
	
	/** Number of sequencing passes performed */
	float passes;
	/** Number of complete sequencing passes */
	int fullPasses=0;
	/** Number of subreads detected in this ZMW */
	int subreads=1;
	/** Count of missing bases or adapters */
	int missing=0;
	/** Number of adapter sequences detected */
	int adapters=0;
	/** Estimated sequencing error rate for this read */
	float errorRate=0;
}