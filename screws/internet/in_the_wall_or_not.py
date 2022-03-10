



import socket
from root.config.main import rAnk, mAster_rank



def whether_in_the_great_fire_wall():
    """Return True if we are in the great cyber fire wall."""
    assert rAnk == mAster_rank, "Should only call it in master core."
    # we use these typical walled websites to check.
    domain_names = ["www.google.com", "www.facebook.com", "www.twitter.com"]
    IN_THE_WALL = [True for _ in range(len(domain_names))]
    for i, dn in enumerate(domain_names):
        try:
            # see if we can resolve the host name -- tells us if there is a DNS listening
            host = socket.gethostbyname(dn)
            # connect to the host -- tells us if the host is actually reachable
            s = socket.create_connection((host, 80), 2)
            s.close()
            IN_THE_WALL[i] = False
        except socket.error:
            IN_THE_WALL[i] = True
    return all(IN_THE_WALL)

