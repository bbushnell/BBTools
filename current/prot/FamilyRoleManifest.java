package prot;

import java.nio.ByteBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;

/**
 * Immutable per-rank role binding for the canonical protein-family assigner.
 *
 * <p>The manifest is deliberately separate from the vector family list: tracked
 * ranks keep their historical vector positions, while appended bait ranks take
 * part in the same shortlist/alignment loop but never become network features.
 * The caller supplies the already-bound expected rep-ID order; this loader
 * requires exact rank/count/order agreement and rejects every inconsistent
 * role/network-feature combination.</p>
 *
 * <p>The legacy four-column schema requires tracked rows to be network
 * features. The metadata-bound schema-7 role artifact also permits tracked
 * rows with {@code network_feature=false}; those rows remain competitive and
 * produce typed masked interceptions. Baits always require false, and all
 * tracked rows precede all bait rows in either schema.</p>
 *
 * @author Yoimiya
 */
public final class FamilyRoleManifest {

	public enum Role{TRACKED, BAIT}

	private FamilyRoleManifest(final String[] repIds_, final Role[] roles_, final boolean[] networkFeatures_,
			final int trackedCount_, final int baitCount_){
		repIds=repIds_; roles=roles_; networkFeatures=networkFeatures_;
		trackedCount=trackedCount_; baitCount=baitCount_;
		assert(repOK()) : "FamilyRoleManifest constructor received inconsistent arrays/counts; loader must validate before publication.";
	}

	/** Loads and validates one role manifest against the binding's exact rep-ID order. */
	public static FamilyRoleManifest load(final String path, final String[] expectedRepIds){
		if(path==null || path.length()==0){throw new IllegalArgumentException("Role-manifest path is blank.");}
		if(expectedRepIds==null || expectedRepIds.length==0){
			throw new IllegalArgumentException("Expected role-manifest rep-ID order is empty.");
		}
		final String[] ids=new String[expectedRepIds.length];
		final Role[] roles=new Role[expectedRepIds.length];
		final boolean[] network=new boolean[expectedRepIds.length];
		final HashSet<String> seen=new HashSet<String>(expectedRepIds.length*2);
		final ByteFile bf=ByteFile.makeByteFile(path,false);
		final CharsetDecoder decoder=StandardCharsets.UTF_8.newDecoder()
			.onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		int row=0, tracked=0, bait=0;
		boolean sawBait=false;
		try{
			final byte[] header=bf.nextLine();
			final String expectedHeader="#rank\trep_id\trole\tnetwork_feature";
			if(header==null || !expectedHeader.equals(new String(header,StandardCharsets.UTF_8))){
				throw new IllegalArgumentException("Unexpected role-manifest header in "+path+"; expected '"+expectedHeader+"'.");
			}
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank role-manifest row at data index "+row+" in "+path+".");}
				if(row>=expectedRepIds.length){
					throw new IllegalArgumentException("Role manifest has more rows than the bound roster (expected "+expectedRepIds.length+"): "+path);
				}
				final String decoded=decode(line,decoder,path,row);
				final String[] fields=decoded.split("\\t",-1);
				if(fields.length!=4){
					throw new IllegalArgumentException("Role-manifest row width is "+fields.length+" instead of 4 at data index "+row+" in "+path+".");
				}
				final int rank=parseRank(fields[0],path,row);
				if(rank!=row){throw new IllegalArgumentException("Role-manifest rank is not contiguous: expected "+row+", got "+rank+" in "+path+".");}
				final String repId=fields[1], roleText=fields[2], networkText=fields[3];
				if(repId.length()==0 || !seen.add(repId)){
					throw new IllegalArgumentException("Blank or duplicate role-manifest rep_id at rank "+rank+" in "+path+": '"+repId+"'.");
				}
				if(!repId.equals(expectedRepIds[rank])){
					throw new IllegalArgumentException("Role-manifest rep_id mismatch at rank "+rank+": expected '"+
						expectedRepIds[rank]+"', got '"+repId+"' in "+path+".");
				}
				final Role role;
				if(roleText.equals("tracked")){role=Role.TRACKED;}
				else if(roleText.equals("bait")){role=Role.BAIT;}
				else{throw new IllegalArgumentException("Unknown role '"+roleText+"' at rank "+rank+" in "+path+".");}
				final boolean feature;
				if(networkText.equals("true")){feature=true;}
				else if(networkText.equals("false")){feature=false;}
				else{throw new IllegalArgumentException("network_feature must be true or false at rank "+rank+" in "+path+".");}
				if((role==Role.TRACKED)!=feature){
					throw new IllegalArgumentException("Role/network_feature mismatch at rank "+rank+" in "+path+
						": role="+roleText+" network_feature="+networkText+".");
				}
				if(role==Role.BAIT){sawBait=true; bait++;}
				else{
					if(sawBait){throw new IllegalArgumentException("Tracked rank follows a bait rank at rank "+rank+" in "+path+".");}
					tracked++;
				}
				ids[rank]=repId; roles[rank]=role; network[rank]=feature; row++;
			}
		}finally{bf.close();}
		if(row!=expectedRepIds.length){
			throw new IllegalArgumentException("Role-manifest row count "+row+" != bound roster count "+expectedRepIds.length+": "+path);
		}
		if(tracked<1){throw new IllegalArgumentException("Role manifest contains no tracked families: "+path);}
		return new FamilyRoleManifest(ids,roles,network,tracked,bait);
	}

	/**
	 * Loads the explicit schema-7 role artifact against the exact family-ID and
	 * representative-ID order of the bound roster.
	 *
	 * <p>Unlike the historical four-column manifest, schema 7 deliberately permits
	 * a tracked family to have {@code network_feature=false}.  Those rows still
	 * participate in competitive assignment but do not occupy neural-network
	 * inputs.  Baits, when present in a later artifact, must also have
	 * {@code network_feature=false}.</p>
	 */
	public static FamilyRoleManifest loadSchema7(final String path,
			final int[] expectedFamilyIds, final String[] expectedRepIds,
			final String expectedRosterSha80, final long expectedMembers){
		if(path==null || path.length()==0){throw new IllegalArgumentException("Role-manifest path is blank.");}
		if(expectedFamilyIds==null || expectedRepIds==null || expectedRepIds.length==0 ||
				expectedFamilyIds.length!=expectedRepIds.length){
			throw new IllegalArgumentException("Schema-7 role binding requires equal nonempty family-ID and rep-ID arrays.");
		}
		DigestSuffix.requireSuffix(expectedRosterSha80,"schema-7 role roster sha80");
		if(expectedMembers<1){throw new IllegalArgumentException("Schema-7 role expected member count must be positive.");}
		final int n=expectedRepIds.length;
		final String[] ids=new String[n];
		final Role[] roles=new Role[n];
		final boolean[] network=new boolean[n];
		final HashSet<Integer> seenFamilyIds=new HashSet<Integer>(n*2);
		final HashSet<String> seenRepIds=new HashSet<String>(n*2);
		final HashMap<String,String> metadata=new HashMap<String,String>();
		final ByteFile bf=ByteFile.makeByteFile(path,false);
		final CharsetDecoder decoder=StandardCharsets.UTF_8.newDecoder()
			.onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		int row=0,tracked=0,bait=0,enabled=0;
		boolean sawColumns=false,sawBait=false;
		final structures.ByteBuilder table=new structures.ByteBuilder(1<<20);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank schema-7 role-manifest line in "+path+".");}
				if(line[0]=='#'){
					if(sawColumns){throw new IllegalArgumentException("Late schema-7 role metadata in "+path+".");}
					final String decoded=decode(line,decoder,path,row);
					final int tab=decoded.indexOf('\t');
					if(tab<2 || tab==decoded.length()-1){throw new IllegalArgumentException("Malformed schema-7 role metadata in "+path+": "+decoded);}
					if(metadata.put(decoded.substring(1,tab),decoded.substring(tab+1))!=null){
						throw new IllegalArgumentException("Duplicate schema-7 role metadata in "+path+": "+decoded.substring(1,tab));
					}
					continue;
				}
				final String decoded=decode(line,decoder,path,row);
				if(!sawColumns){
					if(!"active_index\tfamily_id\trep_id\trole\tnetwork_feature".equals(decoded)){
						throw new IllegalArgumentException("Unexpected schema-7 role column header in "+path+": "+decoded);
					}
					sawColumns=true;
					table.append(line).nl();
					continue;
				}
				table.append(line).nl();
				if(row>=n){throw new IllegalArgumentException("Schema-7 role manifest has more than "+n+" rows: "+path);}
				final String[] fields=decoded.split("\\t",-1);
				if(fields.length!=5){throw new IllegalArgumentException("Schema-7 role row width is "+fields.length+" instead of 5 at rank "+row+" in "+path+".");}
				final int rank=parseRank(fields[0],path,row);
				final int familyId=parseRank(fields[1],path,row);
				if(rank!=row || familyId!=expectedFamilyIds[row] || !seenFamilyIds.add(Integer.valueOf(familyId))){
					throw new IllegalArgumentException("Schema-7 role rank/family mismatch at rank "+row+" in "+path+".");
				}
				final String repId=fields[2];
				if(repId.length()==0 || !repId.equals(expectedRepIds[row]) || !seenRepIds.add(repId)){
					throw new IllegalArgumentException("Schema-7 role rep_id mismatch at rank "+row+" in "+path+".");
				}
				final Role role;
				if("tracked".equals(fields[3])){role=Role.TRACKED; tracked++;}
				else if("bait".equals(fields[3])){role=Role.BAIT; bait++; sawBait=true;}
				else{throw new IllegalArgumentException("Unknown schema-7 role at rank "+row+" in "+path+": "+fields[3]);}
				if(role==Role.TRACKED && sawBait){throw new IllegalArgumentException("Tracked schema-7 role follows bait at rank "+row+" in "+path+".");}
				final boolean feature;
				if("true".equals(fields[4])){feature=true; enabled++;}
				else if("false".equals(fields[4])){feature=false;}
				else{throw new IllegalArgumentException("Schema-7 network_feature must be true or false at rank "+row+" in "+path+".");}
				if(role==Role.BAIT && feature){throw new IllegalArgumentException("Schema-7 bait cannot be a network feature at rank "+row+" in "+path+".");}
				ids[row]=repId; roles[row]=role; network[row]=feature; row++;
			}
		}finally{bf.close();}
		if(!sawColumns || row!=n){throw new IllegalArgumentException("Schema-7 role row count "+row+" != expected "+n+": "+path);}
		requireMetadata(metadata,"artifact_type","schema7_family_roles",path);
		requireMetadata(metadata,"schema","1",path);
		requireMetadata(metadata,"roster_sha80",expectedRosterSha80,path);
		requireMetadata(metadata,"families",Integer.toString(n),path);
		requireMetadata(metadata,"members",Long.toString(expectedMembers),path);
		requireMetadata(metadata,"enabled",Integer.toString(enabled),path);
		requireMetadata(metadata,"masked",Integer.toString(n-enabled),path);
		requireMetadata(metadata,"role_table_sha80_contract",
			"column_header_lf_plus_rows_each_lf",path);
		requireMetadata(metadata,"role_table_sha80",DigestSuffix.bytes(table.toBytes()),path);
		if(metadata.size()!=9){throw new IllegalArgumentException("Schema-7 role metadata count "+metadata.size()+" != 9 in "+path+".");}
		if(tracked<1){throw new IllegalArgumentException("Schema-7 role manifest contains no tracked families: "+path);}
		return new FamilyRoleManifest(ids,roles,network,tracked,bait);
	}

	private static void requireMetadata(final HashMap<String,String> metadata,
			final String key, final String expected, final String path){
		final String observed=metadata.get(key);
		if(!expected.equals(observed)){
			throw new IllegalArgumentException("Schema-7 role metadata #"+key+"='"+observed+"' != '"+expected+"' in "+path+".");
		}
	}

	private static String decode(final byte[] line, final CharsetDecoder decoder, final String path, final int row){
		try{return decoder.decode(ByteBuffer.wrap(line)).toString();}
		catch(CharacterCodingException e){
			throw new IllegalArgumentException("Malformed UTF-8 in role manifest at data index "+row+" in "+path+".",e);
		}
	}

	private static int parseRank(final String text, final String path, final int row){
		if(text.length()==0 || (text.length()>1 && text.charAt(0)=='0')){
			throw new IllegalArgumentException("Noncanonical role-manifest rank '"+text+"' at data index "+row+" in "+path+".");
		}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			if(c<'0' || c>'9'){
				throw new IllegalArgumentException("Noncanonical role-manifest rank '"+text+"' at data index "+row+" in "+path+".");
			}
		}
		try{return Integer.parseInt(text);}
		catch(NumberFormatException e){throw new IllegalArgumentException("Role-manifest rank overflow at data index "+row+" in "+path+".",e);}
	}

	public int size(){return roles.length;}
	public int trackedCount(){return trackedCount;}
	public int baitCount(){return baitCount;}
	public String repId(final int rank){checkRank(rank); return repIds[rank];}
	public Role role(final int rank){checkRank(rank); return roles[rank];}
	public boolean isTracked(final int rank){return role(rank)==Role.TRACKED;}
	public boolean isBait(final int rank){return role(rank)==Role.BAIT;}
	public boolean isNetworkFeature(final int rank){checkRank(rank); return networkFeatures[rank];}

	private void checkRank(final int rank){
		if(rank<0 || rank>=roles.length){throw new IndexOutOfBoundsException("Role-manifest rank "+rank+" outside [0,"+roles.length+").");}
	}

	private boolean repOK(){
		if(repIds==null || roles==null || networkFeatures==null || repIds.length<1 ||
				repIds.length!=roles.length || roles.length!=networkFeatures.length || trackedCount<1 || baitCount<0 ||
				trackedCount+baitCount!=roles.length){return false;}
		boolean sawBait=false;
		for(int i=0; i<roles.length; i++){
			if(repIds[i]==null || repIds[i].length()==0 || roles[i]==null){return false;}
			if(roles[i]==Role.BAIT && networkFeatures[i]){return false;}
			if(roles[i]==Role.BAIT){sawBait=true;}
			else if(sawBait){return false;}
		}
		return true;
	}

	private final String[] repIds;
	private final Role[] roles;
	private final boolean[] networkFeatures;
	private final int trackedCount, baitCount;
}
