package icecream;

import stream.Read;
import structures.ByteBuilder;

/**
 * Mutable sequence and synthetic truth metadata used by IceCreamMaker.
 * Accumulates subread/adapter fragments and emits the tab-separated header consumed
 * by IceCreamGrader. Contamination classification uses the stored subread count.
 * Serialization always uses movie prefix m1_2_3. Each mutable builder belongs to
 * one simulator worker; callers own any sequence arrays shared with it.
 * @author Brian Bushnell, Sucrose
 */
class ReadBuilder{
	
	/**
	 * Wraps the supplied bases without copying them, then initializes metadata.
	 * @param bases_ Non-null sequence array, initially shared with the ByteBuilder
	 * @param passes_ Fractional template passes represented by this fragment
	 * @param movieStart_ Stored starting coordinate
	 * @param zmw_ Stored well identifier
	 */
	public ReadBuilder(byte[] bases_, float passes_, int movieStart_, long zmw_){
		this(new ByteBuilder(bases_), passes_, movieStart_, zmw_);
	}
	
	/**
	 * Borrows the mutable builder and initializes the stored stop to start+length.
	 * A passes value at least one initializes fullPasses to one, not its integer part.
	 * Caller maintains compatible coordinates and metadata; no range checks occur here.
	 * @param bases_ Non-null mutable sequence builder, retained directly
	 * @param passes_ Fractional template passes represented by this fragment
	 * @param movieStart_ Stored starting coordinate
	 * @param zmw_ Stored well identifier
	 */
	public ReadBuilder(ByteBuilder bases_, float passes_, int movieStart_, long zmw_){
		bases=bases_;
		passes=passes_;
		movieStart=movieStart_;
		movieStop=movieStart+bases.length();
		zmw=zmw_;

		fullPasses=passes<1 ? 0 : 1;
	}
	
	/**
	 * Reads the synthetic truth label from the fourth tab-separated header field.
	 * Uses the value after '=' by position without checking the field name; does not
	 * examine sequence bases or accept arbitrary real-read headers.
	 * @param id Non-null synthetic header containing the ordered subreads field
	 * @return Whether the stored subread count is greater than one
	 */
	public static boolean isIceCream(String id){
		String[] terms=id.split("\t");
		int subreads=Integer.parseInt(terms[3].split("=")[1]);
		return subreads>1;
	}
	
	/**
	 * Parses the positional header format emitted by toHeader, retaining r.bases
	 * through a wrapping ByteBuilder. Requires the first six fields; a missing
	 * seventh errorRate field defaults to zero. Field names, coordinate consistency,
	 * and the original movie prefix are not validated or retained separately.
	 * @param r Non-null read with non-null bases and a compatible synthetic header
	 * @return New metadata object whose sequence initially shares the read's array
	 */
	public static ReadBuilder parse(Read r){
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
	
	/** @return The newly formatted synthetic header as a String. */
	@Override
	public String toString(){
		return toHeader().toString();
	}
	
	/**
	 * Serializes m1_2_3/well/start_stop followed by passes, fullPasses, subreads,
	 * missing, adapters and errorRate, in that tab-separated order. Formats passes
	 * with two decimal places and errorRate with three; stored numeric fields are
	 * not recomputed from sequence content.
	 * @return Newly allocated header builder
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
	public int length(){
		return bases.length();
	}
	
	/**
	 * Appends the source sequence, extends stored movieStop by its emitted length,
	 * and sums missing/adapters/fullPasses/subreads/passes. Retains this object's
	 * zmw, movieStart and errorRate; no cross-fragment compatibility check is made.
	 * @param rb Non-null fragment whose metadata the caller considers compatible
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
	 * Creates a Read from a copy of the current bases and the synthetic header.
	 * Passes null qualities and numeric ID zero to the Read constructor.
	 * @return New read with sequence independent of the builder
	 */
	Read toRead(){
		ByteBuilder id=toHeader();
		Read r=new Read(bases.toBytes(), null, id.toString(), 0);
		return r;
	}
	
	/** Mutable sequence builder, borrowed by constructors and mutated by add. */
	ByteBuilder bases;

	/** Stored well identifier, retained when appending fragments. */
	final long zmw;
	/** Stored starting coordinate, retained when appending fragments. */
	final int movieStart;
	/** Stored stop, initialized and extended using emitted sequence lengths. */
	int movieStop;
	
	/** Accumulated fractional template passes, not recomputed from emitted length. */
	float passes;
	/** Complete-pass metadata, initially at most one and summed by add. */
	int fullPasses=0;
	/** Subread fragments represented, initially one; adapter-only fragments set this to zero. */
	int subreads=1;
	/** Missing-adapter joins recorded by the simulator, not a missing-base count. */
	int missing=0;
	/** Included adapter fragments recorded by the simulator, not sequence detections. */
	int adapters=0;
	/** Stored requested simulation error rate; add retains the receiver's value. */
	float errorRate=0;
}
