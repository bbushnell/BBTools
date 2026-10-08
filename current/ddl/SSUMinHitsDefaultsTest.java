package ddl;

/**
 * Guards FindSSU/SSUServer's shared DDL seed-hit default.
 *
 * @author Yelan
 */
public class SSUMinHitsDefaultsTest {

	public static void main(String[] args){
		assert(SSUCompare.DEFAULT_MIN_HITS==5) : SSUCompare.DEFAULT_MIN_HITS;
		assert(new SSUServer(new String[0]).minHits()==5);
		assert(new SSUServer(new String[] {"minhits=8"}).minHits()==8);
		assert(new SSUServer(new String[] {"minhits=5"}).minHits()==5);
		System.err.println("SSUMinHitsDefaultsTest passed.");
	}
}
