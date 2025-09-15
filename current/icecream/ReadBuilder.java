package icecream;

import stream.Read;
import structures.ByteBuilder;

class ReadBuilder {
	
	public ReadBuilder(byte[] bases_, float passes_, int movieStart_, long zmw_) {
		this(new ByteBuilder(bases_), passes_, movieStart_, zmw_);
	}
	
	public ReadBuilder(ByteBuilder bases_, float passes_, int movieStart_, long zmw_) {
		bases=bases_;
		passes=passes_;
		movieStart=movieStart_;
		movieStop=movieStart+bases.length();
		zmw=zmw_;

		fullPasses=passes<1 ? 0 : 1;
	}
	
	public static boolean isIceCream(String id){
		String[] terms=id.split("\t");
		int subreads=Integer.parseInt(terms[3].split("=")[1]);
		return subreads>1;
	}
	
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
	
	public int length() {
		return bases.length();
	}
	
	void add(ReadBuilder rb){
		bases.append(rb.bases);
		
		movieStop+=rb.length();
		missing+=rb.missing;
		adapters+=rb.adapters;
		fullPasses+=rb.fullPasses;
		subreads+=rb.subreads;
		passes+=rb.passes;
	}
	
	Read toRead() {
		//Example: m54283_190403_183820/4194374/919_2614
		//Run ID is m54283_190403_183820
		//zmw ID is 4194374.
		//Read start/stop coordinates are 919_2614
		
		ByteBuilder id=toHeader();
		Read r=new Read(bases.toBytes(), null, id.toString(), 0);
		return r;
	}
	
	ByteBuilder bases;

	final long zmw;
	final int movieStart;
	int movieStop;
	
	float passes;
	int fullPasses=0;
	int subreads=1;
	int missing=0;
	int adapters=0;
	float errorRate=0;
}