package tracker;

import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.Tools;

public class SealStats {
	
//	#File	NYUTS.fq.gz
//	#Total	54848862	8282178162
//	#Matched	936762	1.70790%
//	#Name	Reads	ReadsPct	Bases	BasesPct
//	NYUTS	914178	1.66672%	138040878	1.66672%
	
	public SealStats(String fname_) {
		fname=fname_;
		load(fname);
	}
	
	void load(String fname) {
		ByteFile bf=ByteFile.makeByteFile1(fname, true);
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			lp.set(line);
			if(Tools.startsWith(line, '#')) {//Header
				if(lp.startsWith("#File")) {
					fname=lp.parseString(1);
				}else if(lp.startsWith("#Total")) {
					totalReads=lp.parseLong(1);
					totalBases=lp.parseLong(2);
				}else if(lp.startsWith("#Matched")) {
					matchedReads=lp.parseLong(1);
				}
			}else {
				String name=lp.parseString(0);
				long reads=lp.parseLong(1);
				long bases=lp.parseLong(3);
				SealStatsLine ssl=map.get(name);
				assert(ssl==null);
				ssl=new SealStatsLine(name, reads, bases);
				matchedBases+=bases; //Since that is not reported earlier
				map.put(name, ssl);
			}
		}
		bf.close();
	}
	
	public SealStatsLine countNonmatching(String name) {
		SealStatsLine sum=new SealStatsLine("!"+name, 0, 0);
		for(Entry<String, SealStatsLine> e : map.entrySet()) {
			String key=e.getKey();
			if(!key.equals(name)) {
				SealStatsLine ssl=e.getValue();
				sum.add(ssl);
			}
		}
		return sum;
	}
	
	public SealStatsLine countNonprimary(String name) {
		SealStatsLine primary=primary();
		return countNonmatching(primary.name);
	}
	
	public SealStatsLine primary() {
		SealStatsLine primary=null;
		for(Entry<String, SealStatsLine> e : map.entrySet()) {
			SealStatsLine ssl=e.getValue();
			if(primary==null || ssl.compareTo(primary)>0) {primary=ssl;}
		}
		return primary;
	}
	
	public String fnamePrefix() {
		return ReadWrite.stripToCore(fname);
	}
	
	public static class SealStatsLine implements Comparable<SealStatsLine> {
		
		public SealStatsLine(String name_, long reads_, long bases_) {
			name=name_;
			reads=reads_;
			bases=bases_;
		}
		
		public void add(SealStatsLine ssl) {
			reads+=ssl.reads;
			bases+=ssl.bases;
		}
		
		@Override
		public int compareTo(SealStatsLine o) {
			if(reads!=o.reads) {return reads>o.reads ? 1 : -1;} 
			if(bases!=o.bases) {return bases>o.bases ? 1 : -1;} 
			return name.compareTo(o.name);
		}

		@Override
		public boolean equals(Object o) {return equals((SealStatsLine)o);}
		public boolean equals(SealStatsLine o) {return name.equals(o.name);}
		
		public String name;
		public long reads;
		public long bases;
	}
	
	public String fname;
	public long totalReads, totalBases;
	public long matchedReads, matchedBases;
	public LinkedHashMap<String, SealStatsLine> map=new LinkedHashMap<String, SealStatsLine>();
	
}
