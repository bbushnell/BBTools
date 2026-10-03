package prot;

/**
 * Pluggable ID-hash function for {@link HbmMemberIndexTable} (Increment 3B design v3 sec2 /
 * v5's injectable-hash requirement: the table's correctness must not depend on hash quality
 * -- identity is always decided by exact byte comparison at the probed slot, never the hash
 * value alone. Tests inject a deliberately-colliding hasher against this same interface to
 * prove that property directly, rather than assume it).
 *
 * @author Eru
 */
public interface HbmIdHasher {
	long hash(byte[] idBytes, int off, int len);
}
