import ROOT
import fedrarootlogon
import numpy as np

def np_tracks(track):
    x1 = track.GetSegmentFirst().X()
    y1 = track.GetSegmentFirst().Y()
    z1 = track.GetSegmentFirst().Z()
    x2 = track.GetSegmentLast().X()
    y2 = track.GetSegmentLast().Y()
    z2 = track.GetSegmentLast().Z()
    return(np.array([[x1,y1,z1],[x2,y2,z2]]))

def find_clusters(tracks, proximity_threshold=300.0):
    clusters = []
    visited = set()
    for i, track1 in enumerate(tracks):
        print(i)
        if i in visited:
            continue
        cluster = [np_tracks(track1)]

        visited.add(i)
        for j, track2 in enumerate(tracks):
            if j in visited:
                continue
            dist_start = np.linalg.norm(np_tracks(track1)[0][:2] - np_tracks(track2)[0][:2])
            if dist_start < proximity_threshold:
                cluster.append(np_tracks(track2))
                visited.add(j)
                # print(f"cluster of {len(cluster)}")
                # print(np_tracks(track1))
                # print(np_tracks(track2))
                # print(dist_start)
        if len(cluster) > 2:  # Only keep clusters with more than 2 tracks
            clusters.append(cluster)

    return clusters

def find_staircase_start(cluster):
    # Sort the tracks by the starting z-coordinate
    sorted_tracks = sorted(cluster, key=lambda track: track[0][2])

    direction_vector = sorted_tracks[0][-1] - sorted_tracks[0][0]
    # Check if there is a staircase pattern between any two consecutive tracks
    for i in range(1, len(sorted_tracks)):
        # Check if the current track starts further along the direction vector
        prev_track_start = sorted_tracks[i-1][0]
        curr_track_start = sorted_tracks[i][0]
        if not (np.dot(curr_track_start, direction_vector) > np.dot(prev_track_start, direction_vector)):
            return None, None
    dist_max = np.linalg.norm(sorted_tracks[-1][-1][:2] - sorted_tracks[0][0][:2])
    return sorted_tracks[i-1][0], dist_max  # Return the start position of the first track in the staircase


# Example Usage:

brickid = 121
# f = ROOT.TFile.Open(f"b{brickid:06}.100.0.0.trk.root")
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
scancond.SetSigma0(50,50,0.005,0.005) #change sigma0
scancond.SetDegrad(4) #change angular degradation
gAli.SetScanCond(scancond)

cut = "nseg>4&&t.eScanID.ePlate>5"
dproc.ReadTracksTree(gAli, f"b{brickid:06}.100.0.0.trk.root", cut)
tracks = gAli.eTracks
# [
#     # Each track is represented as a list of two points [start, end]
#     np.array([[0, 0, 0], [1, 1, 1]]),
#     np.array([[2, 2, 2], [3, 3, 3]]),
#     np.array([[1, 1, 0.5], [2, 2, 1.5]]),
#     np.array([[4, 4, 4], [5, 5, 5]]),
# ]

# Find clusters of tracks
clusters = find_clusters(tracks)

# Check each cluster for staircase pattern and print the first track's start position
for cluster in clusters:
    start_position, dist_max = find_staircase_start(cluster)
    if start_position is not None:
        print(f"Staircase found with {len(cluster)} tracks. Maximuim dist = {dist_max}. First track starts at: {start_position}")
