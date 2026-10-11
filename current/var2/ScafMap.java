package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.SamReadInputStream;
import stream.Streamer;
import stream.StreamerFactory;

/**
 * Reference-scaffold registry with insertion-order numeric IDs and name lookup.
 * Loads metadata from SAM/VCF headers and attaches caller-owned reference bases.
 * Exact names take precedence over aliases ending at the first whitespace.
 * Distinct reference names can make an alias ambiguous; header-first promotion
 * has separate reuse rules described by addScaffold(Read).
 * The registry and process-wide default cache are mutable and unsynchronized;
 * finish construction before sharing for lookup and coordinate later mutations.
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Sucrose
 */
public class ScafMap{

	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty registry without changing the process-wide default. */
	public ScafMap(){}

	/*--------------------------------------------------------------*/
	/*----------------      SAM Header Loading      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads scaffold information from a SAM file header.
	 * @param fname SAM filename
	 * @return ScafMap populated with scaffold information
	 */
	@Deprecated
	public static ScafMap loadSamHeader(String fname){return loadSamHeader(fname, null);}

	/**
	 * Loads scaffold information from a SAM file header.
	 * @param ff SAM FileFormat object
	 * @return ScafMap populated with scaffold information
	 */
	@Deprecated
	public static ScafMap loadSamHeader(FileFormat ff){return loadSamHeader(ff, null);}

	/**
	 * Loads scaffold information from a SAM file header into an existing ScafMap.
	 * @param fname SAM filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information
	 */
	@Deprecated
	public static ScafMap loadSamHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, false, false);
		return loadSamHeader(ff, scafMap);
	}

	/**
	 * Starts a one-read SAM/BAM streamer that requests header publication, merges
	 * the JVM-wide shared header, then closes the streamer on normal return.
	 * Shared-header ownership must be coordinated by the caller; this method does
	 * not isolate the header from other inputs or close in a finally block.
	 * @param ff SAM/BAM input descriptor
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information
	 */
	@Deprecated
	public static ScafMap loadSamHeader(FileFormat ff, ScafMap scafMap){
		//maxReads=1, not 0: the streamer's read loop is gated by readID<maxReads and readID only counts
		//READ lines, never @ header lines (SamStreamerST.processFileDirectly). maxReads=0 makes the guard
		//false on entry, so the loop never runs, no @SQ line is ever added, and an EMPTY shared header is
		//published -> callers see scafMap.size()==0 (bbduk vcf=/var= over SAM crashed here). maxReads=1 reads
		//the full header + exactly one read (not the whole file), so @SQ is captured cheaply for SAM and BAM.
		Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, 0, true, false, 1, false);
		ss.start();
		scafMap=waitForSamHeader(scafMap);
		ss.close();
		return scafMap;
	}

	/** Obtains the JVM-wide header with waiting requested and adds its SQ lines.
	 * Does not start a reader. Requires a publisher and nonempty header lines;
	 * the shared-header API may return null when no input has been marked present.
	 * Stops at the first line not beginning with '@'.
	 * @param scafMap Registry to extend, or null to create one
	 * @return Extended or newly created registry */
	public static ScafMap waitForSamHeader(ScafMap scafMap){
		ArrayList<byte[]> header=SamReadInputStream.getSharedHeader(true);
		assert(header!=null);
		if(scafMap==null){scafMap=new ScafMap();}
		for(byte[] line : header){
			if(Tools.startsWith(line, "@SQ\t")){
				scafMap.add(line);
			}else if(line[0]!='@'){
				break;
			}
		}
		return scafMap;
	}

	/*--------------------------------------------------------------*/
	/*----------------      VCF Header Loading      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads scaffold information from a VCF file header.
	 * @param fname VCF filename
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadVcfHeader(String fname){return loadVcfHeader(fname, null);}

	/**
	 * Loads scaffold information from a VCF file header.
	 * @param ff VCF FileFormat object
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadVcfHeader(FileFormat ff){return loadVcfHeader(ff, null);}

	/**
	 * Loads scaffold information from a VCF file header into an existing ScafMap.
	 * @param fname VCF filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadVcfHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.VCF, null, false, false);
		return loadVcfHeader(ff, scafMap);
	}

	/**
	 * Reads leading header lines, adding those beginning with {@code ##contig=<ID=}.
	 * Stops at EOF, an empty line, or the first non-header line. The parser expects
	 * the field layout described by addFromVcf. Closes the ByteFile on normal return;
	 * exceptions propagate without finally-based cleanup.
	 * @param ff VCF FileFormat object
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @return ScafMap populated with scaffold information
	 */
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

	/**
	 * Loads scaffold information from a FASTA reference file.
	 * @param fname FASTA filename
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadReference(String fname, boolean makeDefault){return loadReference(fname, null, null, makeDefault);}

	/**
	 * Loads scaffold information from a FASTA reference file.
	 * @param ff FASTA FileFormat object
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadReference(FileFormat ff, boolean makeDefault){return loadReference(ff, null, null, makeDefault);}

	/**
	 * Loads scaffold information from a FASTA reference file with optional filtering.
	 * @param fname FASTA filename
	 * @param scafMap Existing ScafMap to add to (null to create new)
	 * @param samFilter Filter for scaffold names (null for no filtering)
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information
	 */
	public static ScafMap loadReference(String fname, ScafMap scafMap, SamFilter samFilter, boolean makeDefault){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return loadReference(ff, scafMap, samFilter, makeDefault);
	}

	/**
	 * Reuses the default registry when its stored filename equals ff.name(),
	 * ignoring map, samFilter and makeDefault on that path. This is a filename-only
	 * cache; it does not check whether the input or requested filter has changed.
	 * Otherwise materializes the input reads, filters their IDs, and attaches their
	 * base arrays without copying. Filtering occurs after all reads are loaded.
	 * @param ff FASTA FileFormat object
	 * @param map Existing ScafMap to add to (null to create new)
	 * @param samFilter Filter for scaffold names (null for no filtering)
	 * @param makeDefault Whether to set as the default ScafMap
	 * @return ScafMap populated with scaffold information
	 */
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

	/** Returns the default ScafMap instance.
	 * @return Default ScafMap or null if not set */
	public static ScafMap defaultScafMap(){return defaultScafMap;}

	/**
	 * Sets the default ScafMap instance for global access.
	 * Also disables {@code RNAME_AS_BYTES} in SamLine for string-based lookups.
	 * Does not copy the registry. Asserts that a nonnull filename differs from the
	 * previous nonnull cached filename; publication is not synchronized here.
	 * @param map Registry to store, including null to clear the reference
	 * @param fname Cache key, or null to disable filename-based reuse
	 */
	public static void setDefaultScafMap(ScafMap map, String fname){
		SamLine.RNAME_AS_BYTES=false;
		assert(fname==null || defaultScafMapFile==null || !fname.equals(defaultScafMapFile)) : fname;
		defaultScafMap=map;
		defaultScafMapFile=fname;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Collection Methods     ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes all registry entries and alias bookkeeping. Does not reset the static
	 * default pointer/filename or mutate Scaffold objects held by other callers. */
	public void clear(){
		map.clear();
		alt.clear();
		ambiguousAlt.clear();
		if(promotedAliases!=null){promotedAliases.clear();}
		hasNameOnlyScaffolds=false;
		list.clear();
	}

	/** Returns the number of scaffolds in this ScafMap.
	 * @return Scaffold count */
	public int size(){return list.size();}

	/** Returns the live primary-map key view, which can include promoted full names.
	 * Treat it as read-only: removing keys would leave other registry structures stale.
	 * @return Backed key set, not a snapshot */
	public Set<String> keySet(){
		return map.keySet();
	}

	/** Returns the live view of currently unambiguous first-whitespace aliases.
	 * Treat it as read-only; mutations do not update the primary map or scaffold list.
	 * @return Backed alias key set, not a snapshot */
	public Set<String> altKeySet(){
		return alt.keySet();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Adders             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Attaches the read's base array without copying, replacing any prior bases.
	 * Reuses an exact name, or may promote a prior name-only scaffold through its
	 * first-whitespace alias when trimming is enabled. Once promoted, subsequent
	 * full names sharing that alias can reuse the same object if lengths match;
	 * the stored name and number remain those of the original Scaffold.
	 * Otherwise creates a new scaffold with the next numeric ID.
	 * Length mismatches throw; sequence equality is not checked before replacement.
	 * @param r Nonnull read with a name and nonnull reference bases
	 * @return Added or existing Scaffold object
	 */
	public Scaffold addScaffold(Read r){
		Scaffold scaf=map.get(r.id);
		if(scaf==null && TRIM_WHITESPACE_ALSO && hasNameOnlyScaffolds){
			scaf=findNameOnlyScaffoldForReference(r);
		}
		if(scaf==null){
			scaf=new Scaffold(r.id, size(), r.length());
			add(scaf);
		}
		if(scaf.length!=r.length()){
			throw new IllegalArgumentException("Reference scaffold '"+r.id+"' has length "+r.length()+
					", but the matching scaffold '"+scaf.name+"' has length "+scaf.length+".");
		}
		scaf.bases=r.bases;
		assert(scaf.bases.length==scaf.length) : "Incorrect reference length: "+scaf.bases.length+", "+scaf.length+", "+scaf.name;
		return scaf;
	}

	/**
	 * Parses metadata through Scaffold's fixed-order SQ parser, then resolves an
	 * existing exact name or eligible alias. A reused scaffold keeps its existing
	 * length without comparing the header's length. Otherwise appends metadata only.
	 * @param line Header bytes beginning with @SQ, followed by SN and LN fields
	 * @return Added or existing Scaffold object
	 */
	public Scaffold add(byte[] line){
		Scaffold scaf=new Scaffold(line, size());
		Scaffold old=findScaffoldForSecondarySource(scaf.name);
		if(old!=null){return old;}
		hasNameOnlyScaffolds=true;
		return add(scaf);
	}

	/**
	 * Parses the layout {@code ##contig=<ID=name,length=len>} (also accepts Length).
	 * Takes the name from ID= through the comma preceding length; this is not a
	 * general order-independent VCF attribute parser. Length uses consecutive
	 * decimal digits. Resolves existing exact names or eligible aliases, retaining
	 * the old length without comparing it, or appends a new metadata-only scaffold.
	 * @param line Header bytes with ID immediately followed by length/Length
	 * @return Added or existing Scaffold object
	 */
	public Scaffold addFromVcf(byte[] line){
		String s=new String(line);
		int idStart=s.indexOf("ID=")+3;
		int lengthKeyStart=s.indexOf("length=");
		if(lengthKeyStart<0){lengthKeyStart=s.indexOf("Length=");}
		//Name ends at the comma before "length="
		int commaBeforeLength=s.lastIndexOf(',', lengthKeyStart);
		String name=s.substring(idStart, commaBeforeLength);
		//Parse length from after '=' to first non-digit
		int lengthValStart=s.indexOf('=', lengthKeyStart)+1;
		int lengthValEnd=lengthValStart;
		while(lengthValEnd<line.length && line[lengthValEnd]>='0' && line[lengthValEnd]<='9'){
			lengthValEnd++;
		}
		int length=Parse.parseInt(line, lengthValStart, lengthValEnd);
		Scaffold scaf=new Scaffold(name, size(), length);
		Scaffold old=findScaffoldForSecondarySource(scaf.name);
		if(old!=null){return old;}
		hasNameOnlyScaffolds=true;
		return add(scaf);
	}

	/**
	 * Adds a scaffold with specified name and length.
	 * Reuses only an exact primary-map key, retaining its old length unchecked.
	 * Unlike header adders, this overload does not resolve aliases before adding.
	 * @param s Scaffold name
	 * @param len Scaffold length
	 * @return Added or existing Scaffold object
	 */
	public Scaffold add(String s, int len){
		Scaffold scaf=map.get(s);
		if(scaf!=null){return scaf;}
		scaf=new Scaffold(s, size(), len);
		hasNameOnlyScaffolds=true;
		return add(scaf);
	}

	/**
	 * Appends a scaffold whose number equals the current size and whose exact name
	 * is absent. When trimming is enabled, registers the first-whitespace alias,
	 * including names without whitespace. A collision removes that alias from alt
	 * and marks it ambiguous; exact primary-map entries remain available.
	 * @param scaf Scaffold to add
	 * @return The added scaffold
	 */
	private Scaffold add(Scaffold scaf){
		assert(!map.containsKey(scaf.name));
		assert(size()==scaf.number);

		list.add(scaf);
		map.put(scaf.name, scaf);
		String s=scaf.name;
		if(TRIM_WHITESPACE_ALSO){
			String alias=Tools.trimToWhitespace(s);
			if(!ambiguousAlt.contains(alias)){
				Scaffold old=alt.get(alias);
				if(old!=null && old!=scaf){
					ambiguousAlt.add(alias);
					alt.remove(alias);
				}else{
					alt.put(alias, scaf);
				}
			}
		}
		return scaf;
	}

	/** Resolves a reference through a previously added name-only entry.
	 * First promotion removes the short primary key and inserts the full name;
	 * the Scaffold itself retains its original name/number. Later names with the
	 * promoted alias reuse that object after a length check. Does not compare bases.
	 * @param r Reference read whose full name was not found exactly
	 * @return Reused scaffold, or null if no promotable entry exists */
	private Scaffold findNameOnlyScaffoldForReference(Read r){
		//TODO: Confirm whether promoted aliases may replace different same-length reference sequences; current tests use identical bases (VAR2-006).
		final String fullName=r.id;
		final String alias=Tools.trimToWhitespace(fullName);
		Scaffold scaf=null;

		if(promotedAliases!=null && promotedAliases.contains(alias)){
			scaf=alt.get(alias);
			if(scaf==null){scaf=map.get(alias);}
		}else if(!alias.equals(fullName)){
			Scaffold candidate=map.get(alias);
			if(candidate!=null && candidate.bases==null && candidate.name.equals(alias)){
				if(candidate.length!=r.length()){
					throw lengthCollision(fullName, r.length(), alias, candidate.length);
				}
				if(promotedAliases==null){promotedAliases=new HashSet<String>();}
				promotedAliases.add(alias);
				map.remove(alias);
				scaf=candidate;
			}
		}

		if(scaf!=null){
			if(scaf.length!=r.length()){
				throw lengthCollision(fullName, r.length(), alias, scaf.length);
			}
			map.put(fullName, scaf);
		}
		return scaf;
	}

	/** Describes a full reference name conflicting with existing short-name length. */
	private static IllegalArgumentException lengthCollision(String fullName, int fullLength,
			String shortName, int shortLength){
		return new IllegalArgumentException("Reference scaffold '"+fullName+"' has length "+fullLength+
				", but shortened scaffold name '"+shortName+"' already identifies a scaffold of length "+
				shortLength+". The short name is ambiguous.");
	}

	/** Resolves a SAM/VCF header name by exact primary key, then an exact alias.
	 * Does not trim a longer header name to force a match. A bare ambiguous alias
	 * throws; an unknown name returns null so the caller can add metadata.
	 * @param name Parsed header name
	 * @return Existing scaffold, or null */
	private Scaffold findScaffoldForSecondarySource(String name){
		Scaffold value=map.get(name);
		if(value!=null || !TRIM_WHITESPACE_ALSO){return value;}
		value=alt.get(name);
		if(value!=null){return value;}

		String alias=Tools.trimToWhitespace(name);
		if(!alias.equals(name)){return null;}
		if(ambiguousAlt.contains(alias)){
			throw new IllegalArgumentException("Ambiguous scaffold name '"+alias+
					"' matches multiple reference sequences. Use unique first tokens.");
		}
		return alt.get(alias);
	}

	/**
	 * Adds coverage information from a SAM alignment.
	 * Only processes mapped alignments.
	 * @param sl SAM line containing alignment information
	 */
	public void addCoverage(SamLine sl){
		if(!sl.mapped()){return;}
		Scaffold scaf=getScaffold(sl);
		scaf.add(sl);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Gets scaffold number by name.
	 * @param s Nonnull scaffold name
	 * @return Scaffold number, or -1 for an absent name when assertions are disabled
	 * @throws AssertionError If the name is absent and assertions are enabled
	 * @throws IllegalArgumentException If alias lookup encounters ambiguity
	 */
	public int getNumber(String s){
		Scaffold value=getScaffold(s);
		return value==null ? -1 : value.number;
	}

	/**
	 * Gets scaffold name by number.
	 * @param number Nonnegative scaffold number; negative values throw
	 * @return Scaffold name, or null when number is at least size()
	 */
	public String getName(int number){
		return number>=size() ? null : list.get(number).name;
	}

	/**
	 * Gets scaffold length by number.
	 * @param number Nonnegative scaffold number; negative values throw
	 * @return Scaffold length, or zero when number is at least size()
	 */
	public int getLength(int number){
		return number>=size() ? 0 : list.get(number).length;
	}

	/**
	 * Gets scaffold by numeric ID.
	 * @param number Nonnegative scaffold number; negative values throw
	 * @return Scaffold object, or null when number is at least size()
	 */
	public Scaffold getScaffold(int number){
		return number>=size() ? null : list.get(number);
	}

	/**
	 * Tests whether a scaffold name is present without asserting on absence.
	 * Ambiguous first-token aliases still throw because they identify multiple
	 * possible scaffolds.
	 * @param s Nonnull reference name
	 * @return True when s resolves to a scaffold
	 */
	public boolean containsScaffold(final String s){
		return findScaffold(s, false)!=null;
	}

	/**
	 * Tries an exact primary key, then an already registered exact alias even if
	 * trimming has since been disabled. If still absent and trimming is enabled,
	 * checks alias ambiguity and tries the prefix ending at the first whitespace.
	 * Exact primary keys take precedence over ambiguity in the alias table.
	 * @param s Nonnull reference name
	 * @return Scaffold, or null for an absent name when assertions are disabled
	 * @throws AssertionError If no matching scaffold exists and assertions are enabled
	 * @throws IllegalArgumentException If fallback lookup encounters an ambiguous alias
	 */
	public Scaffold getScaffold(final String s){
		Scaffold value=findScaffold(s, true);
		assert(value!=null) : s+"\n"+keySet()+"\n"+altKeySet()+"\n";
		return value;
	}

	/**
	 * Looks up a scaffold name by exact primary key, exact alias, then optional
	 * first-whitespace alias.
	 * @param s Nonnull reference name
	 * @param assertAbsent Whether to assert if no first-whitespace alias resolves
	 * @return Resolved scaffold, or null when absent
	 */
	private Scaffold findScaffold(final String s, final boolean assertAbsent){
		Scaffold value=map.get(s);
		if(value==null){value=alt.get(s);}
		if(value==null && TRIM_WHITESPACE_ALSO){
			String alias=Tools.trimToWhitespace(s);
			if(ambiguousAlt.contains(alias)){
				throw new IllegalArgumentException("Ambiguous scaffold name '"+alias+
						"' matches multiple reference sequences. Use unique first tokens.");
			}
			if(!alias.equals(s)){value=alt.get(alias);}
			assert(!assertAbsent || value!=null) : "Scaffold not present in reference: "+s+
				"\nalias="+alias+"\nkeySet="+keySet()+"\naltSet="+altKeySet()+"\n";
		}
		return value;
	}

	/**
	 * Gets scaffold from a SAM line’s reference name.
	 * @param sl SAM line
	 * @return Scaffold object
	 */
	public Scaffold getScaffold(SamLine sl){
		return getScaffold(sl.rnameS());
	}

	/**
	 * Gets coverage at a variant position.
	 * @param v Variant object
	 * @return Coverage value
	 */
	public int getCoverage(Var v){
		Scaffold scaf=getScaffold(v.scafnum);
		return scaf.calcCoverage(v);
	}

	/** Calculates total length of all scaffolds.
	 * @return Sum of all scaffold lengths */
	public long lengthSum(){
		long sum=0;
		for(Scaffold scaf : list){
			sum+=scaf.length;
		}
		return sum;
	}

	/** Clears each scaffold's coverage after its readers and writers have finished.
	 * Preserves names, numeric IDs, reference bases and the default-cache binding. */
	public void clearCoverage(){
		for(Scaffold scaf : list){scaf.clearCoverage();}
	}

	/** Returns string representation of all scaffolds.
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

	/** Unique scaffolds in numeric-ID order. */
	final ArrayList<Scaffold> list=new ArrayList<Scaffold>();
	/** Exact lookup keys; several promoted full names may reference one Scaffold. */
	final HashMap<String, Scaffold> map=new HashMap<String, Scaffold>();
	/** First-whitespace aliases that have not collided with a different Scaffold. */
	private final HashMap<String, Scaffold> alt=new HashMap<String, Scaffold>();
	/** Colliding aliases removed from alt; reset only by clear(). */
	private final HashSet<String> ambiguousAlt=new HashSet<String>();
	/** Short metadata names promoted by reference loading; lazily allocated. */
	private HashSet<String> promotedAliases;
	/** Records that a metadata-only adder inserted an entry since the last clear. */
	private boolean hasNameOnlyScaffolds=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared default reference registry, not copied or synchronized. */
	private static ScafMap defaultScafMap=null;
	/** Exact filename cache key; null disables loadReference's cache hit. */
	private static String defaultScafMapFile=null;
	/** Controls alias creation and fallback lookup; changing it does not rebuild aliases. */
	public static boolean TRIM_WHITESPACE_ALSO=true;

	/** Legacy identifier; this class does not implement Serializable. */
	private static final long serialVersionUID=1L;
}
