package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Parse;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;

/**
 * Maps scaffold (chromosome/contig) names to Scaffold objects and provides
 * efficient lookup by name or numeric ID. Supports loading scaffold information
 * from SAM headers, VCF headers, or FASTA reference files.
 * 
 * Maintains both primary and alternative name mappings to handle whitespace
 * variations in scaffold names across different file formats.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class ScafMap {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an empty ScafMap */
	public ScafMap(){}

	/*--------------------------------------------------------------*/
	/*----------------      SAM Header Loading      ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads scaffold information from SAM file header
	 * @param fname SAM filename
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadSamHeader(String fname){return loadSamHeader(fname, null);}
	
	/** Loads scaffold information from SAM file header
	 * @param ff SAM FileFormat object
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadSamHeader(FileFormat ff){return loadSamHeader(ff, null);}
	
	/** Loads scaffold information from SAM file header into existing ScafMap
	 * @param fname SAM filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadSamHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, false, false);
		return loadSamHeader(ff, scafMap);
	}
	
	/** Loads scaffold information from SAM file header into existing ScafMap
	 * @param ff SAM FileFormat object
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadSamHeader(FileFormat ff, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(ff);
		if(scafMap==null){scafMap=new ScafMap();}
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			if(Tools.startsWith(line, "@SQ\t")){
				scafMap.add(line);
			}else if(line[0]!='@'){
				break;
			}
			line=bf.nextLine();
		}
		bf.close();

		return scafMap;
	}

	/*--------------------------------------------------------------*/
	/*----------------      VCF Header Loading      ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads scaffold information from VCF file header
	 * @param fname VCF filename
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadVcfHeader(String fname){return loadVcfHeader(fname, null);}
	
	/** Loads scaffold information from VCF file header
	 * @param ff VCF FileFormat object
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadVcfHeader(FileFormat ff){return loadVcfHeader(ff, null);}
	
	/** Loads scaffold information from VCF file header into existing ScafMap
	 * @param fname VCF filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadVcfHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.VCF, null, false, false);
		return loadVcfHeader(ff, scafMap);
	}
	
	/** Loads scaffold information from VCF file header into existing ScafMap
	 * @param ff VCF FileFormat object
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadVcfHeader(FileFormat ff, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(ff);
		if(scafMap==null){scafMap=new ScafMap();}
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			if(Tools.startsWith(line, "##contig=<ID=")){
				scafMap.addFromVcf(line);
			}else if(line[0]!='#'){
				break;
			}
			line=bf.nextLine();
		}
		bf.close();
		return scafMap;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Reference File Loading    ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads scaffold information from FASTA reference file
	 * @param fname FASTA filename
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadReference(String fname, boolean makeDefault){return loadReference(fname, null, null, makeDefault);}
	
	/** Loads scaffold information from FASTA reference file
	 * @param ff FASTA FileFormat object
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadReference(FileFormat ff, boolean makeDefault){return loadReference(ff, null, null, makeDefault);}
	
	/** Loads scaffold information from FASTA reference file with filtering
	 * @param fname FASTA filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @param samFilter Filter for scaffold names (null for no filtering)
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadReference(String fname, ScafMap scafMap, SamFilter samFilter, boolean makeDefault){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return loadReference(ff, scafMap, samFilter, makeDefault);
	}
	
	/** Loads scaffold information from FASTA reference file with filtering
	 * @param ff FASTA FileFormat object
	 * @param map Existing ScafMap to add to (null to create new)
	 * @param samFilter Filter for scaffold names (null for no filtering)
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information */
	public static ScafMap loadReference(FileFormat ff, ScafMap map, SamFilter samFilter, boolean makeDefault){
		if(defaultScafMapFile!=null && defaultScafMapFile.equals(ff.name())){return defaultScafMap;}
		ArrayList<Read> reads=ReadInputStream.toReads(ff, -1);
		if(map==null){map=new ScafMap();}
		for(Read r : reads){
			if(samFilter==null || samFilter.passesFilter(r.id)){
				map.addScaffold(r);
			}
		}
		if(makeDefault){setDefaultScafMap(map, ff.name());}
		return map;
	}

	/** Returns the default ScafMap instance
	 * @return Default ScafMap or null if not set */
	public static ScafMap defaultScafMap(){return defaultScafMap;}
	
	/** Sets the default ScafMap instance for global access
	 * @param map ScafMap to set as default
	 * @param fname Filename associated with this ScafMap */
	public static void setDefaultScafMap(ScafMap map, String fname){
		SamLine.RNAME_AS_BYTES=false;
		assert(fname==null || defaultScafMapFile==null || !fname.equals(defaultScafMapFile)) : fname;
		defaultScafMap=map;
		defaultScafMapFile=fname;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Collection Methods     ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Removes all scaffolds from this ScafMap */
	public void clear(){
		map.clear();
		alt.clear();
		list.clear();
	}
	
	/** Returns the number of scaffolds in this ScafMap
	 * @return Scaffold count */
	public int size(){return list.size();}
	
	/** Returns set of primary scaffold names
	 * @return Set of scaffold names */
	public Set<String> keySet() {
		return map.keySet();
	}
	
	/** Returns set of alternative scaffold names
	 * @return Set of alternative names */
	public Set<String> altKeySet() {
		return alt.keySet();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Adders             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Adds a scaffold from a Read object, storing sequence data
	 * @param r Read containing scaffold name and sequence
	 * @return Added or existing Scaffold object */
	public Scaffold addScaffold(Read r){
		Scaffold scaf=map.get(r.id);
		if(scaf==null){
			scaf=new Scaffold(r.id, size(), r.length());
			add(scaf);
		}
		scaf.bases=r.bases;
		assert(scaf.bases.length==scaf.length) : "Incorrect reference length: "+scaf.bases.length+", "+scaf.length+", "+scaf.name;
		return scaf;
	}
	
	/** Adds a scaffold from SAM header line (@SQ format)
	 * @param line SAM header line
	 * @return Added or existing Scaffold object */
	public Scaffold add(byte[] line){
		Scaffold scaf=new Scaffold(line, size());
		Scaffold old=map.get(scaf.name);
		if(old!=null){return old;}
		return add(scaf);
	}
	
	/** Adds a scaffold from VCF contig header line (##contig=<ID=name,length=len>)
	 * @param line VCF contig header line
	 * @return Added or existing Scaffold object */
	public Scaffold addFromVcf(byte[] line){
		int comma=Tools.indexOf(line, (byte)',');
		String name=new String(line, 13, comma-13);
		int stop=Tools.indexOf(line, (byte)',', comma+1);
		stop=(stop<0 ? line.length : stop);
		int length=Parse.parseInt(line, comma+8, stop-1);
		Scaffold scaf=new Scaffold(name, size(), length);
		Scaffold old=map.get(scaf.name);
		if(old!=null){return old;}
		return add(scaf);
	}
	
	/** Adds a scaffold with specified name and length
	 * @param s Scaffold name
	 * @param len Scaffold length
	 * @return Added or existing Scaffold object */
	public Scaffold add(String s, int len){
		Scaffold scaf=map.get(s);
		if(scaf!=null){return scaf;}
		scaf=new Scaffold(s, size(), len);
		return add(scaf);
	}
	
	/** Internal method to add a scaffold and update all lookup structures
	 * @param scaf Scaffold to add
	 * @return The added scaffold */
	private Scaffold add(Scaffold scaf){
		assert(!map.containsKey(scaf.name));
		assert(size()==scaf.number);
		
		list.add(scaf);
		map.put(scaf.name, scaf);
		String s=scaf.name;
		if(TRIM_WHITESPACE_ALSO){
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					String s2=s.substring(0, i);
					boolean b=alt.containsKey(s2);
					assert(!b);
					if(!b){
						alt.put(s2, scaf);
					}
					break;
				}
			}
		}
		return scaf;
	}
	
	/** Adds coverage information from a SAM alignment
	 * @param sl SAM line containing alignment information */
	public void addCoverage(SamLine sl){
		if(!sl.mapped()){return;}
		Scaffold scaf=getScaffold(sl);
		scaf.add(sl);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets scaffold number by name
	 * @param s Scaffold name
	 * @return Scaffold number or -1 if not found */
	public int getNumber(String s){
		Scaffold value=getScaffold(s);
		return value==null ? -1 : value.number;
	}
	
	/** Gets scaffold name by number
	 * @param number Scaffold number
	 * @return Scaffold name or null if invalid number */
	public String getName(int number){
		return number>=size() ? null : list.get(number).name;
	}
	
	/** Gets scaffold length by number
	 * @param number Scaffold number
	 * @return Scaffold length or 0 if invalid number */
	public int getLength(int number){
		return number>=size() ? 0 : list.get(number).length;
	}
	
	/** Gets scaffold by numeric ID
	 * @param number Scaffold number
	 * @return Scaffold object or null if invalid number */
	public Scaffold getScaffold(int number){
		return number>=size() ? null : list.get(number);
	}
	
	/** Gets scaffold by name with intelligent whitespace handling
	 * @param s Scaffold name (may contain whitespace variations)
	 * @return Scaffold object */
	public Scaffold getScaffold(final String s){
		Scaffold value=map.get(s);
		if(value==null){value=alt.get(s);}
		if(value==null && TRIM_WHITESPACE_ALSO){
			int index=-1;
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					index=i;
					break;
				}
			}
			String sub=null;
			if(index>0){
				sub=s.substring(0, index);
				value=alt.get(sub);
				if(value==null){value=map.get(sub);}
			}
			assert(value!=null) : "Scaffold not present in reference: "+s+"\nwhitespace="+index+
				"\nsubstring="+sub+"\nkeySet="+keySet()+"\naltSet="+altKeySet()+"\n";
		}
		assert(value!=null) : s+"\n"+keySet()+"\n"+altKeySet()+"\n";
		return value;
	}
	
	/** Gets scaffold from SAM line reference name
	 * @param sl SAM line
	 * @return Scaffold object */
	public Scaffold getScaffold(SamLine sl){
		return getScaffold(sl.rnameS());
	}
	
	/** Gets coverage at a variant position
	 * @param v Variant object
	 * @return Coverage value */
	public int getCoverage(Var v){
		Scaffold scaf=getScaffold(v.scafnum);
		return scaf.calcCoverage(v);
	}
	
	/** Calculates total length of all scaffolds
	 * @return Sum of all scaffold lengths */
	public long lengthSum() {
		long sum=0;
		for(Scaffold scaf : list){
			sum+=scaf.length;
		}
		return sum;
	}
	
	/** Clears coverage information from all scaffolds */
	public void clearCoverage() {
		for(Scaffold scaf : list){scaf.clearCoverage();}
	}
	
	/** Returns string representation of all scaffolds
	 * @return Formatted scaffold information */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		for(Scaffold sc : list){sb.append(sc).append('\n');}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Ordered list of scaffolds (indexed by scaffold number) */
	final ArrayList<Scaffold> list=new ArrayList<Scaffold>();
	/** Primary name-to-scaffold mapping */
	final HashMap<String, Scaffold> map=new HashMap<String, Scaffold>();
	/** Alternative name-to-scaffold mapping (for whitespace variants) */
	private final HashMap<String, Scaffold> alt=new HashMap<String, Scaffold>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global default ScafMap instance */
	private static ScafMap defaultScafMap=null;
	/** Filename associated with default ScafMap */
	private static String defaultScafMapFile=null;
	/** Whether to create alternative mappings for whitespace-trimmed names */
	public static boolean TRIM_WHITESPACE_ALSO=true;
	
	/** Serialization ID */
	private static final long serialVersionUID = 1L;
}