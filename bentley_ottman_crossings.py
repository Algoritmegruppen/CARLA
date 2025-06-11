import heapq


class Event:
    def __init__(self, x, segments, is_start):
        self.x = x
        self.segments = segments
        self.is_start = is_start

    def __lt__(self, other):
        return self.x < other.x


class Segment:
    def __init__(self, start, end):
        if start[0] > end[0]:
            start, end = end, start
        self.start = start
        self.end = end

    def __lt__(self, other):
        return self.start < other.start

    def intersects(self, other):
        def ccw(a, b, c):
            return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

        return ccw(self.start, other.start, other.end) != ccw(
            self.end, other.start, other.end
        ) and ccw(self.start, self.end, other.start) != ccw(
            self.start, self.end, other.end
        )


def bentley_ottmann(G, pos):
    # Convert edges of the graph to segments
    segments = [Segment(pos[u], pos[v]) for u, v in G.edges()]

    # Initialize event queue
    events = []
    for seg in segments:
        heapq.heappush(events, Event(seg.start[0], [seg], True))
        heapq.heappush(events, Event(seg.end[0], [seg], False))

    active_segments = []
    crossing_count = 0

    while events:
        event = heapq.heappop(events)

        if event.is_start:
            for seg in active_segments:
                for new_seg in event.segments:
                    if (
                        tuple(seg.start) == tuple(new_seg.start)
                        or tuple(seg.start) == tuple(new_seg.end)
                        or tuple(seg.end) == tuple(new_seg.end)
                        or tuple(seg.end) == tuple(new_seg.start)
                    ):
                        continue
                    if seg.intersects(new_seg):
                        crossing_count += 1

            active_segments.extend(event.segments)
        else:
            for seg in event.segments:
                try:
                    active_segments.remove(seg)
                except:
                    pass

    return crossing_count
