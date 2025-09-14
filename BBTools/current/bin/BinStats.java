package bin;

public class BinStats implements Comparable<BinStats> {
	
	BinStats(Bin b, String name_){
		name=name_;
		id=b.id();
		taxid=b.taxid;
		if(taxid<1) {taxid=b.labelTaxid;}
		size=b.size();
		contigs=b.numContigs();
		contam=b.contam;
		complt=b.completeness;
		badContigs=b.badContigs;
		gc=b.gc();
		depth=b.depth();
		minDepth=b.minContigDepth();
		maxDepth=b.maxContigDepth();

		assert(b.gc()!=0) : this;
	}
	
	@Override
	public int compareTo(BinStats b) {
		if(size!=b.size) {return size>b.size ? -1 : 1;}
		return name.compareTo(b.name);
	}
	
	String type(boolean useRNA) {
		return type(complt, contam, r16Scount, r23Scount, r5Scount, trnaCount, useRNA);
	}

	public boolean hq(boolean useRNA) {
		String type=type(useRNA);
		return type.endsWith("HQ");
	}
	public boolean mq(boolean useRNA) {
		String type=type(useRNA);
		return type.equals("MQ");
	}
	
//	static String type(float complt, float contam) {
//		if(contam<0.01 && complt>=0.99) {return "UHQ";}
//		if(contam<0.02 && complt>=0.95) {return "VHQ";}
//		if(contam<0.05 && complt>=0.90) {return "HQ";}
//		if(contam<0.10 && complt>=0.50) {return "MQ";}
//		if(contam<0.20 || complt<0.20) {return "VLQ";}
//		return "LQ";
//	}
	
	static String type(float complt, float contam, int r16S, int r23S, int r5S, int trna, boolean useRNA) {
		boolean rnaOK=!useRNA || (r16S>0 && r23S>0 && r5S>0 && trna>=18);
		if(rnaOK) {
			if(contam<0.01 && complt>=0.99) {return "UHQ";}
			if(contam<0.02 && complt>=0.95) {return "VHQ";}
			if(contam<0.05 && complt>=0.90) {return "HQ";}
		}
		if(contam<0.10 && complt>=0.50) {return "MQ";}
		if(contam<0.20 || complt<0.20) {return "VLQ";}
		return "LQ";
	}
	
	final String name;
	int id;
	int taxid;
	long size;
	int contigs;
	int badContigs;
	float contam;
	float complt;
	float gc;
	float depth;
	float minDepth, maxDepth;

	int r5Scount=0;
	int r16Scount=0;
	int r18Scount=0;
	int r23Scount=0;
	int trnaCount=0;
	int cdsCount=0;
	long cdsLength=0;
	public String lineage;
	
}
